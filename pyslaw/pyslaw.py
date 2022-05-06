"""
Library for solving linear atmospheric gravity wave problems
"""
import numpy as np

#TODO: create class with objects. For instance, you could have
# different functions for calculating eta, w, p, u. These functions
# could take z as input (to easily calculate output at 1 height)
# Also, think about whether the solver should be able to handle 2D and 3D
# or if you want separate solvers for that

class LinearModel(object):
    def __init__(self,x,h,U,N):
        # assert x, h are one-dimensional
        # assert x is even
        assert(x.size==h.size)
        assert(np.unique(np.diff(x)).size==1)
        assert(all(np.isreal(h)))

        self.U = U
        self.N = N
        self.Nx = x.size
        self.dx = np.unique(np.diff(x))

        self.k = 2.0 * np.pi * np.fft.rfftfreq(self.Nx,self.dx)
        self.m = self.vertical_wavenumbers()

        self.hc = np.fft.rfft(h)
        

    def calculate_eta(self,z):
        if np.isscalar(z): z = np.array([z])

        etac = self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        # Set defunct modes to zero
        etac[-1,:] = 0.
        return np.squeeze(np.fft.irfft(etac,axis=0))

    def calculate_w(self,z):
        if np.isscalar(z): z = np.array([z])

        wc = 1j*self.U*self.k[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        # Set defunct modes to zero
        wc[-1,:] = 0.
        return np.squeeze(np.fft.irfft(wc,axis=0))

    def calculate_u(self,z):
        if np.isscalar(z): z = np.array([z])

        uc = -1j*self.U*self.m[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        # Set defunct modes to zero
        uc[-1,:] = 0.
        return np.squeeze(np.fft.irfft(uc,axis=0))

    def calculate_p(self,z):
        if np.isscalar(z): z = np.array([z])

        pc = 1j * self.U**2 * self.m[:,np.newaxis] * self.hc[:,np.newaxis] * np.exp(1j*self.m[:,np.newaxis]*z)
        # Set defunct modes to zero
        pc[-1,:] = 0.
        return np.squeeze(np.fft.irfft(pc,axis=0))

    def vertical_wavenumbers(self):
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Propagating waves
        iprop = np.where((self.U*self.k)**2>self.N**2)
        #Evanescent waves (excluding where U*k=0, for which m is set to zero) 
        ievan = np.where(~(((-self.U*self.k)==0) | ((self.U*self.k)**2>self.N**2)))
    
        m[iprop] = 1j*np.abs(self.k[iprop])*np.sqrt(1-self.N**2/(self.U*self.k[iprop])**2)
        m[ievan] = -np.sign(-self.U*self.k[ievan])*np.abs(self.k[ievan])*np.sqrt(self.N**2/(self.U*self.k[ievan])**2-1)
        return m

