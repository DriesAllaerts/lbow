"""
Library for solving transient linear buoyancy wave problems
"""
import numpy as np

class OneLayerModel(object):
    def __init__(self,x,t,h,U,N):
        """Initialize model and set governing parameters

        Args:
            x (array): x coordinates (2D grid)
            t (array): time coordinates (2D grid)
            h (array): surface elevation
            U (float): wind speed
            N (float): Brunt-Vaisala frequency
        """
        # assert x, h are one-dimensional
        #assert(x.size % 2 == 0), 'Size of x must be even'
        #assert(x.size == h.size), 'Size of x must match size of h'
        #assert(np.unique(np.diff(x)).size==1), 'x must be spaced equidistantly'
        #assert(all(np.isreal(h))), 'h should be real-valued'
        assert(U != 0), 'Background wind speed should be non-zero'

        # Store wind speed and Brunt Vaisala frequency
        self.U = U
        self.N = N
       
        # Determine if input grid follows 'ij' or 'xy' indexing 
        if np.unique(np.diff(x,axis=0))[0] != 0.:
            indexing = 'ij'
            xaxis = 0
            taxis = 1
        else:
            indexing = 'xy'
            xaxis = 1
            taxis = 0

        # Calculate horizontal wave numbers and frequencies
        dx = np.unique(np.diff(x,axis=xaxis))
        dt = np.unique(np.diff(t,axis=taxis))

        self.Nx = x.shape[xaxis]
        self.Nt = t.shape[taxis]

        # We use rfft2 for efficiency, which means the last axis
        # will be transformed using rfft and accordingly only uses
        # half of the wavenumbers/frequencies
        if indexing == 'ij':
            # last axis corresponds to t
            k = 2.0 * np.pi * np.fft.fftfreq(self.Nx,dx)
            omega = - 2.0 * np.pi * np.fft.rfftfreq(self.Nt,dt)
        else:
            # last axis corresponds to x
            k = 2.0 * np.pi * np.fft.rfftfreq(self.Nx,dx)
            omega = - 2.0 * np.pi * np.fft.fftfreq(self.Nt,dt)
        self.k, self.omega = np.meshgrid(k,omega,indexing=indexing)
        # Note the minus sign for the frequencies. This is because in linear theory
        # the solution is assumed to be a plane wave of the form exp[i(k*x-omega*t)],
        # whereas the 2D Fourier transform uses exp[i(k*x+omega*t)]

        #Intrinsic frequency
        self.Omega = self.omega - self.U*self.k

        # Calculate vertical wave numbers
        self.m = self.vertical_wavenumbers()

        # Store FFT of input signal h
        self.hc = np.fft.rfft2(h,norm='forward')

    def vertical_wavenumbers(self):
        """Calculate vertical wavenumbers

        Returns:
            array: vertical wave numbers
        """
        m = np.zeros(self.k.shape,dtype=np.complex128)
        #Evanescent waves
        ievan = np.where(self.Omega**2>self.N**2)
        #Propagating waves (excluding where Omega=0, for which m is set to zero) 
        iprop = np.where(~((self.Omega==0) | (self.Omega**2>self.N**2)))
    
        m[ievan] = 1j*np.abs(self.k[ievan])*np.sqrt(1-self.N**2/self.Omega[ievan]**2)
        # w_g = -Omega*m/kappa**2, so choose sign(m)=-sign(Omega).
        m[iprop] = -np.sign(self.Omega[iprop])*np.abs(self.k[iprop])*np.sqrt(self.N**2/self.Omega[iprop]**2-1)
        return m
        

class HalfPlaneModel(OneLayerModel):
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
        assert(all(z>=0)), 'All z must be positive'

        if varname == 'eta':
            # Solving d**2 eta/dz**2 + m**2 z = 0
            # Solution of the form eta = A exp(jmz) + B exp(-jmz)
            A = self.hc
        elif varname == 'w':
            # From definition w = D(eta)/Dt
            A = -1j * self.Omega * self.hc
        elif varname == 'u':
            # From continuity equation du/dx + dw/dz = 0
            #A = -1j * self.Omega * self.m * self.hc / self.k
            pass
        elif varname == 'p':
            # From x-momentum equation U * du/dx = - dp/dx
            #A = 1j * self.U * self.Omega * self.m * self.hc / self.k
            pass

        var = A[np.newaxis,...] * np.exp(1j*self.m[np.newaxis,...]*z[:,np.newaxis,np.newaxis])
        # Set defunct modes to zero
#        var[-1,:] = 0.
        return np.squeeze(np.fft.irfft2(var,axes=[-2,-1],norm='forward'))
