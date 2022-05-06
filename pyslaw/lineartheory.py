"""
Library for solving linear atmospheric gravity wave problems
"""
import numpy as np

class LinearModel(object):
    def __init__(self,x,h,U,N):
        # assert x, h are one-dimensional
        assert(x.size % 2 == 0), 'Size of x must be even'
        assert(x.size == h.size), 'Size of x must match size of h'
        assert(np.unique(np.diff(x)).size==1), 'x must be spaced equidistantly'
        assert(all(np.isreal(h))), 'h should be real-valued'

        # Store wind speed and Brunt Vaisala frequency
        self.U = U
        self.N = N
        
        # Calculate horizontal wave numbers
        dx = np.unique(np.diff(x))
        self.k = 2.0 * np.pi * np.fft.rfftfreq(x.size,dx)
        # Calculate vertical wave numbers
        self.m = self.vertical_wavenumbers()

        # Store FFT of input signal h
        self.hc = np.fft.rfft(h)
        

    def solve(self,varname,z):
        assert(varname in ['eta','u','w','p'])

        if np.isscalar(z): z = np.array([z])

        if varname == 'eta':
            var = self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        elif varname == 'w':
            # From definition w = U * d(eta)/dx
            var = 1j*self.U*self.k[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        elif varname == 'u':
            # From continuity equation du/dx + dw/dz = 0
            var = -1j*self.U*self.m[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        elif varname == 'p':
            # From x-momentum equation U * du/dx = - dp/dx
            var = 1j * self.U**2 * self.m[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)

        # Set defunct modes to zero
        var[-1,:] = 0.
        return np.squeeze(np.fft.irfft(var,axis=0))

    def vertical_wavenumbers(self):
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Propagating waves
        iprop = np.where((self.U*self.k)**2>self.N**2)
        #Evanescent waves (excluding where U*k=0, for which m is set to zero) 
        ievan = np.where(~(((-self.U*self.k)==0) | ((self.U*self.k)**2>self.N**2)))
    
        m[iprop] = 1j*np.abs(self.k[iprop])*np.sqrt(1-self.N**2/(self.U*self.k[iprop])**2)
        m[ievan] = -np.sign(-self.U*self.k[ievan])*np.abs(self.k[ievan])*np.sqrt(self.N**2/(self.U*self.k[ievan])**2-1)
        return m

