"""
Library for solving linear atmospheric gravity wave problems
"""
import numpy as np

#TODO: Decide on sign of m in the solution, not when calculating m in general
#TODO: Make iprop and ievan separate functions
class LinearModel(object):
    def __init__(self,x,h,U,N):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates
            h (array): surface elevation at x coordinates
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
        """
        # assert x, h are one-dimensional
        assert(x.size % 2 == 0), 'Size of x must be even'
        assert(x.size == h.size), 'Size of x must match size of h'
        assert(np.unique(np.diff(x)).size==1), 'x must be spaced equidistantly'
        assert(all(np.isreal(h))), 'h should be real-valued'
        assert(U != 0), 'Background wind speed should be non-zero'

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
        """Solve model at specified heights

        Args:
            varname (str): name of the variable to be calculated (eta, u, w, or p)
            z (float or array): heights at which the variable is calculated

        Returns:
            array: model solution at x coordinates and specified heights
        """
        assert(varname in ['eta','u','w','p'])

        if np.isscalar(z): z = np.array([z])

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
        return np.squeeze(np.fft.irfft(var,axis=0))

    def vertical_wavenumbers(self):
        """Calculate vertical wavenumbers

        Returns:
            array: vertical wave numbers
        """
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Evanescent waves
        ievan = np.where((self.U*self.k)**2>self.N**2)
        #Propagating waves (excluding where U*k=0, for which m is set to zero) 
        iprop = np.where(~(((-self.U*self.k)==0) | ((self.U*self.k)**2>self.N**2)))
    
        m[ievan] = 1j*np.abs(self.k[ievan])*np.sqrt(1-self.N**2/(self.U*self.k[ievan])**2)
        # w_g = -Omega*m/kappa**2, so choose sign(m)=-sign(Omega).
        # for stationary waves, omega=Omega+U*k=0, so Omega=-U*k
        m[iprop] = -np.sign(-self.U*self.k[iprop])*np.abs(self.k[iprop])*np.sqrt(self.N**2/(self.U*self.k[iprop])**2-1)
        return m

