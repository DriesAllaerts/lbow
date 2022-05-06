"""
Library for solving linear atmospheric gravity wave problems
"""
import numpy as np

#TODO: create class with objects. For instance, you could have
# different functions for calculating eta, w, p, u. These functions
# could take z as input (to easily calculate output at 1 height)
# Also, think about whether the solver should be able to handle 2D and 3D
# or if you want separate solvers for that

def calculate_agw(x,z,h,U,N):
    # assert x, h are one-dimensional
    # assert x is even
    assert(x.size==h.size)
    assert(np.unique(np.diff(x)).size==1)
    assert(all(np.isreal(h)))
    Nx = x.size
    dx = np.unique(np.diff(x))

    k = 2.0 * np.pi * np.fft.rfftfreq(Nx,dx)

    hc = np.fft.rfft(h)
    m = vertical_wavenumber(k,U,N)

    etac = hc[:,np.newaxis] * np.exp(1j*m[:,np.newaxis]*z)
    # Set defunct modes to zero
    etac[:,-1] = 0.
    return np.fft.irfft(etac,axis=0)

def vertical_wavenumber(k,U,N):
    m = np.zeros(k.shape,dtype=np.complex128)
    #Propagating waves
    iprop = np.where((U*k)**2>N**2)
    #Evanescent waves (excluding where U*k=0, for which m is set to zero) 
    ievan = np.where(~(((-U*k)==0) | ((U*k)**2>N**2)))

    m[iprop] = 1j*np.abs(k[iprop])*np.sqrt(1-N**2/(U*k[iprop])**2)
    m[ievan] = -np.sign(-U*k[ievan])*np.abs(k[ievan])*np.sqrt(N**2/(U*k[ievan])**2-1)
    return m
