(sec:HalfPlaneModel)=
# Half Plane Model

In the Half Plane Model, the flow is solved on an infinite half-plane above the ground, and the upper boundary condition effectively applies at infinity. Different boundary conditions apply for propagating and evanescent waves.

We first derive the solution for the general, transient flow scenario. Next, we present the simpler solution for a steady state flow scenario.

## Transient solution

As waves are only excited at the surface, the solution can only contain waves that propagate/radiate upwards. Following the sign convention for the vertical wave number introduced in {numref}`Section %s <sec:verticalWaveNumber>`, we need to retain the term with $m_1=m$ in the general solution {eq}`eqn:generalSolution` and discard the term with $m_2=-m$. This gives

$$
\hat{\eta}_1(k,z,\omega) = A(k,\omega)e^{jm(k,\omega)z},
$$ (eqn:eta_HalfPlaneModel)

with $A(k,\omega)=\hat{h}(k,\omega)$ following from equation {eq}`eqn:Fourier_BC`, and the vertical wave number for 1D models ($l=0$) given by

$$
    m(k,\omega) = \begin{cases}
    j|k|\sqrt{1-N^2/\Omega^2} & \text{for}\;\Omega^2>N^2 \\
    -\text{sign}(\Omega)\,|k|\sqrt{N^2/\Omega^2-1} & \text{for}\;\Omega^2<N^2
    \end{cases}
$$ (eqn:m_HalfPlaneModel)

With the vertical displacement given by equations {eq}`eqn:eta_HalfPlaneModel` and {eq}`eqn:m_HalfPlaneModel`, the other flow variables can be calculated using the {ref}`polarization equations <sec:polarizationEqns>`:

$$
    \hat{w}_1(k,z,\omega)=-j\Omega\hat{\eta}_1(k,z,\omega)
$$ (eqn:w_HalfPlaneModel)

$$
    \hat{u}_1(k,z,\omega)=\frac{jm\Omega}{k}\hat{\eta}_1(k,z,\omega)
$$ (eqn:u_HalfPlaneModel)

$$
    \hat{p}_1(k,z,\omega)=\rho_0\frac{jm\Omega^2}{k^2}\hat{\eta}_1(k,z,\omega)
$$ (eqn:p_HalfPlaneModel)

The full flow solution for any flow variable is found by applying the 2D inverse Fourier transform:

$$
    \phi_1(x,z,t)=\iint_{-\infty}^{\infty}\hat{\phi}_1(k,z,\omega)\exp{[j(kx-\omega t)]}\;\mathrm{d} k \,\mathrm{d} \omega
$$

with $\phi_1$ representing any variable in $(\eta_1,u_1,w_1,p_1)$.
```{note}
The solution is developed for a plane wave of the form $\sim\exp{[j(kx-\omega t)]}$ whereas the 2D Fourier transform is normally defined using $\sim\exp{[j(kx+\omega t)]}$. This needs to be taken into account in the numerical implementation.
```

## Steady state solution
For steady state flow solutions, equations {eq}`eqn:eta_HalfPlaneModel` to {eq}`eqn:p_HalfPlaneModel` can be further simplified since $\omega=0$ and hence $\Omega=-u_0k$. For completeness, the simplified solution is given below:

$$
\hat{\eta}_1(k,z) = \hat{h}(k)e^{jm(k)z},
$$

with

$$
    \hat{h}(k)=\frac{1}{2\pi}\int_{-\infty}^{\infty}h(x)\exp{(-jkx)}\;\mathrm{d} x
$$ (eqn:h_k)
and

$$
    m(k) = \begin{cases}
    j|k|\sqrt{1-N^2/(u_0^2k^2)} & \text{for}\;u_0^2k^2>N^2, \\
    -\text{sign}(-u_0k)\,|k|\sqrt{N^2/(u_0^2k^2)-1} & \text{for}\;u_0^2k^2<N^2,
    \end{cases}
$$ (eqn:m_k)

Other flow variables are given by

$$
\hat{w}_1(k,z)=j\, u_0 k\,\hat{\eta}_1(k,z),
$$

$$
\hat{u}_1(k,z)=-j\, u_0 m\,\hat{\eta}_1(k,z),
$$

$$
 \hat{p}_1(k,z)=j\, \rho_0 u_0^2 m\,\hat{\eta}_1(k,z).
$$

The full flow solution is given by the 1D inverse Fourier transform

$$
\phi_1(x,z)=\int_{-\infty}^{\infty}\hat{\phi}_1(k,z)\exp{(jkx)}\;\mathrm{d} k
$$

with $\phi_1$ representing any variable in $(\eta_1,u_1,w_1,p_1)$