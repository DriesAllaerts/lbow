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

class OneLayerModel(object):
    """
    Base class for steady state models consisting of one layer
    """
    def __init__(self,x,h,U,N):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates
            h (array): surface elevation at x coordinates
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
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

        # Store FFT of input signal h
        self.hc = np.fft.rfft(h,norm='forward')

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
    def __init__(self,x,h,U,N,H):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates
            h (array): surface elevation at x coordinates
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
            H (float): Height of the top boundary
        """
        super().__init__(x,h,U,N)

        assert(H>0), 'H must be positive'
        self.H = H
    
    def solve(self,varname,z):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated

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

        return np.squeeze(np.fft.irfft(var,axis=0,norm='forward'))

class HalfPlaneModel(OneLayerModel):
    """
    Steady state one-layer model with the ground surface as bottom boundary
    and the top boundary at infinity (radiation boundary condition)
    """
    def solve(self,varname,z):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated

        Returns:
            array: model solution at x coordinates and specified heights
        """
        assert(varname in ['eta','u','w','p'])

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
        return np.squeeze(np.fft.irfft(var,axis=0,norm='forward'))


class MultiLayerModel(object):
    """
    Class for steady state models consisting of multiple layers
    """
    def __init__(self,x,hs,U,N,hi):
        """Initialize model and set governing parameters

        Args:
            x  (array): x coordinates
            hs (array): surface elevation at x coordinates
            U  (array): array of wind speeds in the various layers
            N  (array): array of Brunt-Vaisala frequency in the various layers
            hi (array): heights marking start and end of the various layers
        """

        # Assertions
        assert(len(x.shape)==1), 'x must be a one-dimensional array'
        assert(x.size % 2 == 0), 'size of x must be even'
        assert(x.shape == hs.shape), 'x and hs must have same dimensions'
        assert(all(np.isreal(hs))), 'hs should be real-valued'
        assert(len(U) == len(N)), 'U and N should have the same size'
        assert(len(U) == len(hi)), 'U and hi should have the same size'
        assert(all(np.array(U) != 0)), 'background wind speed should be non-zero'
        assert(all(np.array(N) >= 0)), 'Brunt-Vaisala frequency should be non-negative'
        assert(all(np.diff(hi) > 0)), 'hi should be monotonically increasing'

        dx = np.unique(np.diff(x))
        assert(np.allclose(dx,dx[0])), 'x must be spaced equidistantly'

        # Store wind speed and Brunt Vaisala frequency
        self.Us = U
        self.Ns = N

        # Store heights and number of layers
        self.Nl = len(U)
        self.hi = hi
        
        # Calculate horizontal wave numbers
        self.Nx = x.size
        self.k = 2.0 * np.pi * np.fft.rfftfreq(self.Nx,dx[0])
        self.Nk = self.k.size

        # Calculate vertical wave numbers
        self.m = np.zeros((self.Nk,self.Nl),dtype=np.complex128)
        for l in range(self.Nl):
            self.m[:,l] = self.vertical_wavenumbers(self.Us[l],self.Ns[l])

        # Store FFT of input signal h
        self.hc = np.fft.rfft(hs,norm='forward')

    def vertical_wavenumbers(self,U,N):
        """Calculate vertical wavenumbers

        Returns:
            array: vertical wave numbers
        """
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Evanescent waves
        ievan = np.where((-U*self.k)**2>N**2)
        #Propagating waves (excluding where U*k=0, for which m is set to zero) 
        iprop = np.where(~(((-U*self.k)==0) | ((-U*self.k)**2>N**2)))
    
        m[ievan] = 1j*np.abs(self.k[ievan])*np.sqrt(1-N**2/(-U*self.k[ievan])**2)
        # w_g = -Omega*m/kappa**2, so choose sign(m)=-sign(Omega).
        # for stationary waves, omega=Omega+U*k=0, so Omega=-U*k
        m[iprop] = -np.sign(-U*self.k[iprop])*np.abs(self.k[iprop])*np.sqrt(N**2/(-U*self.k[iprop])**2-1)
        return m

    def solve(self,varname,z):
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): height(s) at which the variable is calculated

        Returns:
            array: model solution at x coordinates and specified heights
        """
        assert(varname in ['eta','u','w','p'])

        if varname in ['u','p']:
            print('Solution in terms of u or p not yet implemented')
            return

        if np.isscalar(z): z = np.array([z])
        assert(all(z>=0)), 'All z must be positive'

        # Create array of heights including a "top" height
        hi = np.zeros(self.Nl+1)
        hi[:-1] = self.hi[:]
        hi[-1] = max(np.max(z)+1.e-6,self.hi[-1])

        # For every k, solve system of linear equations for coefficients
        # --------------------------------------------------------------
        # Initialize A and B matrix
        A = np.zeros((self.Nk,2*self.Nl,2*self.Nl),dtype=np.complex128)
        B = np.zeros((self.Nk,2*self.Nl),dtype=np.complex128)

        # Surface boundary condition
        A[:,0,0] = 1.0
        A[:,0,1] = np.exp(1j*self.m[:,0]*(hi[1]-hi[0]))
        B[:,0] = self.hc

        # Interface boundary conditions
        for l in range(1,self.Nl):
            # Eta continuous
            A[:,2*l-1,2*(l-1)]   =  np.exp(1j*self.m[:,l-1]*(hi[l]-hi[l-1]))
            A[:,2*l-1,2*(l-1)+1] =  1.0
            A[:,2*l-1,2*(l-1)+2] = -1.0
            A[:,2*l-1,2*(l-1)+3] = -np.exp(1j*self.m[:,l]*(hi[l+1]-hi[l]))
            # d(eta)/dz continuous
            A[:,2*l,2*(l-1)]   =  1j*self.m[:,l-1] * np.exp(1j*self.m[:,l-1]*(hi[l]-hi[l-1]))
            A[:,2*l,2*(l-1)+1] = -1j*self.m[:,l-1]
            A[:,2*l,2*(l-1)+2] = -1j*self.m[:,l]
            A[:,2*l,2*(l-1)+3] =  1j*self.m[:,l]   * np.exp(1j*self.m[:,l]*(hi[l+1]-hi[l]))
            # Exception for k=0: Set B_(l-1) to zero
            A[0,2*l,2*(l-1)+1] =  1

        # Top boundary condition
        A[:,-1,-1] = 1

        # Solve linear system for every k
        X = np.linalg.solve(A,B)

        # Compose solution in frequency space
        # ----------------------------------

        var = np.zeros((self.Nk,z.size),dtype=np.complex128)

        if varname == 'eta':
            # eta = A exp( jm (z-zlow) ) + B exp( jm (zhigh-z) )
            for l in range(self.Nl):
                mask = (z >= hi[l] ) & (z < hi[l+1])
                # Ignore overflow, this should be captured by the mask
                with np.errstate(over='ignore'):
                    var += X[:,2*l,np.newaxis]   * np.where(mask[np.newaxis,:],np.exp(1j*self.m[:,l,np.newaxis]*(z-hi[l])),0)
                    var += X[:,2*l+1,np.newaxis] * np.where(mask[np.newaxis,:],np.exp(1j*self.m[:,l,np.newaxis]*(hi[l+1]-z)),0)

        elif varname == 'w':
            # From definition w = U * d(eta)/dx
            for l in range(self.Nl):
                mask = (z >= hi[l] ) & (z < hi[l+1])
                # Ignore overflow, this should be captured by the mask
                with np.errstate(over='ignore'):
                    var += 1j * self.Us[l] * self.k[:, np.newaxis] * X[:,2*l,np.newaxis]   * \
                                np.where(mask[np.newaxis,:],np.exp(1j*self.m[:,l,np.newaxis]*(z-hi[l])),0)
                    var += 1j * self.Us[l] * self.k[:, np.newaxis] * X[:,2*l+1,np.newaxis] * \
                                np.where(mask[np.newaxis,:],np.exp(1j*self.m[:,l,np.newaxis]*(hi[l+1]-z)),0)

        assert(not np.isnan(var).any()), 'Solution contains nans, probably due to overflow'
        # Set defunct modes to zero
        var[-1,:] = 0.
        return np.squeeze(np.fft.irfft(var,axis=0,norm='forward'))



# Build A matrix

#        elif varname == 'u':
#            # From continuity equation du/dx + dw/dz = 0
#            A = -1j*self.U*self.m * self.hc
#        elif varname == 'p':
#            # From x-momentum equation U * du/dx = - dp/dx
#            A = 1j * self.U**2 * self.m * self.hc
