"""
Library for solving two-dimensional steady state linear buoyancy wave problems

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
    Base class for 2D steady state models consisting of one layer
    """
    def __init__(self,x,y,h,U,V,N,hydrostatic=False,fftw_flag='FFTW_ESTIMATE'):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates (2D grid)
            y (array): y coordinates (2D grid)
            h (array): surface elevation
            U (float): wind speed (x-component)
            V (float): wind speed (y-component)
            N (float): buoyancy frequency
            hydrostatic (bool): Make hydrostatic assumption (default to False)
            fftw_flag(string, optional): flag for the fftw algorithm
        """
        # Assertions
        assert(len(x.shape)==2), 'x must be a two-dimensional grid'
        assert(x.shape == y.shape), 'x and y must have same dimensions'
        assert(x.shape == h.shape), 'x and h must have same dimensions'
        assert(np.isreal(h).all()), 'h should be real-valued'
        assert((U^2+V^2) != 0), 'Background wind speed should be non-zero'

        # Store wind speed and buoyancy frequency
        self.U = U
        self.V = V
        self.N = N
       
        # Determine if input grid follows 'ij' or 'xy' indexing 
        if np.unique(np.diff(x,axis=0))[0] != 0.:
            indexing = 'ij'
            xaxis = 0
            yaxis = 1
        else:
            indexing = 'xy'
            xaxis = 1
            yaxis = 0

        # Assert input coordinates have an even number of grid points
        assert(x.shape[xaxis] % 2 == 0), 'Number of grid points in x direction should be even'
        assert(x.shape[yaxis] % 2 == 0), 'Number of grid points in y direction should be even'

        self.Nx = x.shape[xaxis]
        self.Ny = y.shape[yaxis]

        # Calculate horizontal wave number
        dx = np.unique(np.diff(x,axis=xaxis))
        dy = np.unique(np.diff(y,axis=yaxis))

        assert(np.allclose(dx,dx[0])), 'x must be spaced equidistantly'
        assert(np.allclose(dy,dy[0])), 'y must be spaced equidistantly'

        # We use rfft2 for efficiency, which means the last axis
        # will be transformed using rfft and accordingly only uses
        # half of the wavenumbers/frequencies
        if indexing == 'ij':
            # last axis corresponds to y
            k = 2.0 * np.pi * np.fft.fftfreq(self.Nx,dx[0])
            l = 2.0 * np.pi * np.fft.rfftfreq(self.Ny,dy[0])
            self.kDefunctModeIndex = (int(self.Nx/2),slice(None))
            self.lDefunctModeIndex = (slice(None),-1)
        else:
            # last axis corresponds to x
            k = 2.0 * np.pi * np.fft.rfftfreq(self.Nx,dx[0])
            l = 2.0 * np.pi * np.fft.fftfreq(self.Ny,dy[0])
            self.kDefunctModeIndex = (slice(None),-1)
            self.lDefunctModeIndex = (int(self.Ny/2),slice(None))
        self.k, self.l = np.meshgrid(k,l,indexing=indexing)

        #Intrinsic frequency
        self.Omega = - self.U*self.k - self.V*self.l
        self.zeroOmegaIndex = np.isclose(self.Omega/np.max(np.abs(self.Omega)),0,atol=1e-06)

        # Calculate vertical wave numbers
        self.m = self.vertical_wavenumbers(hydrostatic)

        # Set up forward fft routine
        hr = pyfftw.empty_aligned(h.shape,dtype='float64')
        hc = pyfftw.empty_aligned(self.k.shape,dtype='complex128')

        fft_object = pyfftw.FFTW(hr,hc,axes=(0,1),
            flags=(fftw_flag,),
            direction='FFTW_FORWARD',
            threads=multiprocessing.cpu_count(),
            normalise_idft=False
        )

        # Store FFT of input signal h
        self.hc = fft_object(h)

    def vertical_wavenumbers(self,hydrostatic=False):
        """Calculate vertical wavenumbers

        Args:
            hydrostatic (bool): Make hydrostatic assumption (default to False)
        Returns:
            array: vertical wave numbers
        """
        m = np.zeros(self.k.shape,dtype=np.complex128)

        if hydrostatic:
            with np.errstate(divide='ignore',invalid='ignore'):
                m = -np.sqrt(self.k**2+self.l**2)*self.N/self.Omega
            m[self.zeroOmegaIndex] = 0.
        else:
            #Evanescent waves
            ievan = self.Omega**2>self.N**2
            #Propagating waves (excluding where Omega=0, for which m is set to zero) 
            iprop = np.logical_and(~(self.Omega**2>self.N**2), ~self.zeroOmegaIndex)
        
            m[ievan] = 1j*np.sqrt(self.k[ievan]**2+self.l[ievan]**2) * np.sqrt(1-self.N**2/self.Omega[ievan]**2)
            # w_g = -Omega*m/kappa**2, so choose sign(m)=-sign(Omega).
            m[iprop] = -np.sign(self.Omega[iprop])*np.sqrt(self.k[iprop]**2+self.l[iprop]**2) * np.sqrt(self.N**2/self.Omega[iprop]**2-1)
        return m
        

class HalfPlaneModel(OneLayerModel):
    """
    Two-dimensional, steady state one-layer model with the ground surface as bottom boundary
    and the top boundary at infinity (radiation boundary condition)
    """
    def solve(self,varname,z,fftw_flag='FFTW_ESTIMATE'):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated
            fftw_flag(string, optional): flag for the fftw algorithm

        Returns:
            array: model solution at x,y coordinates and specified heights
        """
        assert(varname in ['eta','u','v','w','p'])

        if np.isscalar(z): z = np.array([z])
        assert(all(z>=0)), 'All z must be positive'

        if varname == 'eta':
            # Solving d**2 eta/dz**2 + m**2 z = 0
            # Solution of the form eta = A exp(jmz) + B exp(-jmz)
            A = self.hc
        elif varname == 'w':
            # From definition w = D(eta)/Dt
            A = -1j * self.Omega * self.hc
        elif varname == 'u':
            print('Currently, you can only solve for eta or w, solving for u is not (yet) implemented')
            pass
#            with np.errstate(divide='ignore',invalid='ignore'):
#                A = 1j * self.Omega * self.m * self.hc / self.k
#            A[self.k == 0] = 0.
        elif varname == 'v':
            print('Currently, you can only solve for eta or w, solving for v is not (yet) implemented')
            pass
#            with np.errstate(divide='ignore',invalid='ignore'):
#                A = 1j * self.Omega * self.m * self.hc / self.k
#            A[self.k == 0] = 0.
        elif varname == 'p':
            print('Currently, you can only solve for eta or w, solving for p is not (yet) implemented')
            pass

        var = A[np.newaxis,...] * np.exp(1j*self.m[np.newaxis,...]*z[:,np.newaxis,np.newaxis])

        # Set defunct modes to zero
        var[(slice(None),) + self.kDefunctModeIndex] = 0.0
        var[(slice(None),) + self.lDefunctModeIndex] = 0.0

        # Set up inverse fft routine
        var_c = pyfftw.empty_aligned(var.shape,dtype='complex128')
        var_r = pyfftw.empty_aligned((var.shape[0],var.shape[1],2*var.shape[2]-2),dtype='float64')

        ifft_object = pyfftw.FFTW(var_c,var_r,axes=(-2,-1),
            flags=(fftw_flag,),
            direction='FFTW_BACKWARD',
            threads=multiprocessing.cpu_count(),
            normalise_idft=False
        )

        return np.squeeze(ifft_object(var))
        #return np.squeeze(np.fft.irfft2(var,axes=[-2,-1],norm='forward'))
