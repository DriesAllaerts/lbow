"""
Library for solving steady state linear buoyancy wave problems

Copyright 2022 Dries Allaerts

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np
import multiprocessing
import pyfftw

class OneLayerModel(object):
    """
    Base class for steady state models consisting of one layer
    """
    def __init__(self,x,h,U,N,fftw_flag='FFTW_ESTIMATE'):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates
            h (array): surface elevation at x coordinates
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
            fftw_flag(string, optional): flag for the fftw algorithm
        """

        # Assertions
        assert(len(x.shape)==1), 'x must be a one-dimensional array'
        assert(x.size % 2 == 0), 'size of x must be even'
        assert(x.shape == h.shape), 'x and h must have same dimensions'
        assert(all(np.isreal(h))), 'h should be real-valued'
        assert(U != 0), 'background wind speed should be non-zero'

        dx = np.unique(np.diff(x))
        assert(np.allclose(dx,dx[0])), 'x must be spaced equidistantly'

        # Store wind speed and Brunt Vaisala frequency
        self.U = U
        self.N = N
        
        # Calculate horizontal wave numbers
        self.Nx = x.size
        self.k = 2.0 * np.pi * np.fft.rfftfreq(self.Nx,dx[0])
        # Calculate vertical wave numbers
        self.m = self.vertical_wavenumbers()

        # Set up forward fft routine
        hr = pyfftw.empty_aligned(h.shape,dtype='float64')
        hc = pyfftw.empty_aligned(self.k.shape,dtype='complex128')

        fft_object = pyfftw.FFTW(hr,hc,axes=(0,),
            flags=(fftw_flag,),
            direction='FFTW_FORWARD',
            threads=multiprocessing.cpu_count(),
            normalise_idft=False
        )

        # Store FFT of input signal h
        self.hc = fft_object(h)

    def vertical_wavenumbers(self):
        """Calculate vertical wavenumbers

        Returns:
            array: vertical wave numbers
        """
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Evanescent waves
        ievan = np.where((-self.U*self.k)**2>self.N**2)
        #Propagating waves (excluding where U*k=0, for which m is set to zero) 
        iprop = np.where(~(((-self.U*self.k)==0) | ((-self.U*self.k)**2>self.N**2)))
    
        m[ievan] = 1j*np.abs(self.k[ievan])*np.sqrt(1-self.N**2/(-self.U*self.k[ievan])**2)
        # w_g = -Omega*m/kappa**2, so choose sign(m)=-sign(Omega).
        # for stationary waves, omega=Omega+U*k=0, so Omega=-U*k
        m[iprop] = -np.sign(-self.U*self.k[iprop])*np.abs(self.k[iprop])*np.sqrt(self.N**2/(-self.U*self.k[iprop])**2-1)
        return m
        
class ChannelModel(OneLayerModel):
    """
    Steady state one-layer model with the ground surface as bottom boundary
    and a rigid lid as the top boundary
    """
    def __init__(self,x,h,U,N,H,fftw_flag='FFTW_ESTIMATE'):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates
            h (array): surface elevation at x coordinates
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
            H (float): Height of the top boundary
            fftw_flag(string, optional): flag for the fftw algorithm
        """
        super().__init__(x,h,U,N,fftw_flag)

        assert(H>0), 'H must be positive'
        self.H = H
    
    def solve(self,varname,z,fftw_flag='FFTW_ESTIMATE'):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated
            fftw_flag(string, optional): flag for the fftw algorithm

        Returns:
            array: model solution at x coordinates and specified heights
        """
        assert(varname in ['eta','u','w','p'])

        if np.isscalar(z): z = np.array([z])
        assert(all((z>=0) & (z<=self.H))), 'All z must lie between 0 and H = {:f}'.format(self.H)

        if varname in ['eta', 'w']:
            # Solving d**2 eta/dz**2 + m**2 z = 0
            # Solution of the form eta = A exp(jmz) + B exp(-jmz)

            # For k=0, m=0 which results in division by zero.
            # Ignore division by zero and set mean mode to zero later on
            with np.errstate(divide='ignore',invalid='ignore'):
                eta1 = self.hc[:,np.newaxis] \
                        / (1 - np.exp(1j*2*self.m[:,np.newaxis]*self.H)) \
                        * np.exp(1j*self.m[:,np.newaxis]*z)
                eta2 = -self.hc[:,np.newaxis] \
                        / (1 - np.exp(1j*2*self.m[:,np.newaxis]*self.H)) \
                        * np.exp(1j*self.m[:,np.newaxis]*(2*self.H-z))            
            var = eta1 + eta2


            # From definition w = U * d(eta)/dx
            if varname == 'w': var *= self.U * 1j * self.k[:, np.newaxis]

        elif varname in ['u', 'p']:
            # From continuity equation du/dx + dw/dz = 0

            # For k=0, m=0 which results in division by zero.
            # Ignore division by zero and set mean mode to zero later on
            with np.errstate(divide='ignore',invalid='ignore'):
                u1 = -1j * self.m[:, np.newaxis] * self.U * self.hc[:,np.newaxis] \
                        / (1 - np.exp(1j*2*self.m[:,np.newaxis]*self.H)) \
                        * np.exp(1j*self.m[:,np.newaxis]*z)
                u2 = -1j * self.m[:, np.newaxis] * self.U * self.hc[:,np.newaxis] \
                        / (1 - np.exp(1j*2*self.m[:,np.newaxis]*self.H)) \
                        * np.exp(1j*self.m[:,np.newaxis]*(2*self.H-z))
            var = u1 + u2

            # From x-momentum equation U * du/dx = - dp/dx
            if varname == 'p': var *= -self.U

        # Set mean and defunct modes to zero
        var[-1,:] = 0.
        var[0,:]  = 0.

        # Set up inverse fft routine
        var_c = pyfftw.empty_aligned(var.shape,dtype='complex128')
        var_r = pyfftw.empty_aligned((2*var.shape[0]-2,var.shape[1]),dtype='float64')

        ifft_object = pyfftw.FFTW(var_c,var_r,axes=(0,),
            flags=(fftw_flag,),
            direction='FFTW_BACKWARD',
            threads=multiprocessing.cpu_count(),
            normalise_idft=False
        )

        return np.squeeze(ifft_object(var))
        #return np.squeeze(np.fft.irfft(var,axis=0,norm='forward'))

class HalfPlaneModel(OneLayerModel):
    """
    Steady state one-layer model with the ground surface as bottom boundary
    and the top boundary at infinity (radiation boundary condition)
    """
    def solve(self,varname,z,space='real',fftw_flag='FFTW_ESTIMATE'):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated
            space (str): return solution in 'real' or 'fourier' space
            fftw_flag(string, optional): flag for the fftw algorithm

        Returns:
            array: model solution at x (or k) coordinates and specified heights
        """
        assert(varname in ['eta','u','w','p'])
        assert(space in ['real','fourier'])

        if np.isscalar(z): z = np.array([z])
        assert(all(z>=0)), 'All z must be positive'

        if varname == 'eta':
            # Solving d**2 eta/dz**2 + m**2 z = 0
            # Solution of the form eta = A exp(jmz) + B exp(-jmz)
            A = self.hc
        elif varname == 'w':
            # From definition w = U * d(eta)/dx
            A = 1j*self.U*self.k * self.hc
        elif varname == 'u':
            # From continuity equation du/dx + dw/dz = 0
            A = -1j*self.U*self.m * self.hc
        elif varname == 'p':
            # From x-momentum equation U * du/dx = - dp/dx
            A = 1j * self.U**2 * self.m * self.hc

        var = A[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        # Set defunct modes to zero
        var[-1,:] = 0.

        if space == 'fourier':
            return np.squeeze(var)
        else:
            # Set up inverse fft routine
            var_c = pyfftw.empty_aligned(var.shape,dtype='complex128')
            var_r = pyfftw.empty_aligned((2*var.shape[0]-2,var.shape[1]),dtype='float64')

            ifft_object = pyfftw.FFTW(var_c,var_r,axes=(0,),
                flags=(fftw_flag,),
                direction='FFTW_BACKWARD',
                threads=multiprocessing.cpu_count(),
                normalise_idft=False
            )

            return np.squeeze(ifft_object(var))
            #return np.squeeze(np.fft.irfft(var,axis=0,norm='forward'))
