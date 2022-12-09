(sec:HalfPlaneModel)=
# Half Plane Model

In the Half Plane Model, the flow is solved on an infinite half-plane above the ground, and the upper boundary condition effectively applies at infinity. Different boundary conditions apply for propagating and evanescent waves.

We first derive the solution for the general, transient flow scenario. Next, we present the simpler solution for a steady state flow scenario.

## Transient solution

For evanescent waves ($\Omega^2>N^2$), the negative imaginary root

$$
    m_2=-j|k|\sqrt{1-N^2/\Omega^2}
$$
will give rise to unphysical growth of the perturbation amplitude with height, so the coefficient $B$ needs to be set to zero in order for the solution to remain bounded.

For propagating waves ($\Omega^2<N^2$), a radiation condition applies, which states that as waves are only excited at the surface, the vertical group velocity of the waves at infinity should be directed upward. The group velocity relative to the flow is defined as $\mathbf{c}_g=(\partial\Omega/\partial k,\partial\Omega/\partial l,\partial\Omega/\partial m)$. The group velocity relative to the ground can be found by taking the derivative of the apparent frequency $\omega$ rather than the intrinsic frequency, but since $\omega=\Omega+\mathbf{u}_0\cdot\mathbf{k}$, this comes down to a simple vector summation. As the background vertical velocity is assumed to be zero, $w_g=\partial\Omega/\partial m$. The partial derivative of $\Omega$ to $m$ can be found by reformulating the expression of the vertical wavenumber $m$ (eq.~\ref{eqn:dispersionEquation}) in terms of $\Omega$, which gives $\Omega=\pm N|k|/\sqrt{k^2+m^2}$. Then, it can be shown that $w_g=-\Omega m/(k^2+m^2)$. Hence, for $w_g$ to be positive, the vertical wave number must be chosen so that $\text{sign}(m)=-\text{sign}(\Omega)$.

For both propagating and evanescent waves, the boundary condition implies that only one root of equation {eq}`eqn:dispersionEquation` should be used. Hence, the solution can be summarised as follows:
```{math}
:label: eqn:eta_HalfPlaneModel
\hat{\eta}_1(k,z,\omega) = A(k,\omega)e^{jm(k,\omega)z},
```
with $A(k,\omega)=\hat{h}(k,\omega)$ following from equation {eq}`eqn:Fourier_BC`, and the vertical wave number given by
```{math}
:label: eqn:m_HalfPlaneModel
    m(k,\omega) = \begin{cases}
    j|k|\sqrt{1-N^2/\Omega^2} & \text{for}\;\Omega^2>N^2 \\
    -\text{sign}(\Omega)\,|k|\sqrt{N^2/\Omega^2-1} & \text{for}\;\Omega^2<N^2
    \end{cases}
```

With the vertical displacement given by equations {eq}`eqn:eta_HalfPlaneModel` and {eq}`eqn:m_HalfPlaneModel`, other flow variables like $u_1$, $w_1$, and $p_1$ are easily related to the vertical displacement $\hat{\eta}_1(k,z,\omega)$. The vertical velocity $w_1$ follows from the kinematic condition {eq}`eqn:kinematicCondition`:
```{math}
:label: eqn:w_HalfPlaneModel
    \hat{w}_1(k,z,\omega)=-j\Omega\hat{\eta}_1(k,z,\omega)
```
The streamwise velocity component $u_1$ is obtained from the continuity equation {eq}`eqn:cont`:
```{math}
:label: eqn:u_HalfPlaneModel
    \hat{u}_1(k,z,\omega)=-\frac{m}{k}\hat{w}_1(k,z,\omega)=\frac{jm\Omega}{k}\hat{\eta}_1(k,z,\omega)
```
The pressure follows from the momentum equation {eq}`eqn:umom`:
```{math}
:label: eqn:p_HalfPlaneModel
    \hat{p}_1(k,z,\omega)=\rho_0\frac{\Omega}{k}\hat{u}_1(k,z,\omega)=\rho_0\frac{jm\Omega^2}{k^2}\hat{\eta}_1(k,z,\omega)
```
Finally, the full flow solution for any flow variable is found by applying the 2D inverse Fourier transform:
```{math}
    \phi_1(x,z,t)=\iint_{-\infty}^{\infty}\hat{\phi}_1(k,z,\omega)\exp{[j(kx-\omega t)]}\;\mathrm{d} k \,\mathrm{d} \omega
```
with $\phi_1$ representing any variable in $(\eta_1,u_1,w_1,p_1)$.
```{note}
The solution is developed for a plane wave of the form $\sim\exp{[j(kx-\omega t)]}$ whereas the 2D Fourier transform is normally defined using $\sim\exp{[j(kx+\omega t)]}$. This needs to be taken into account in the numerical implementation.
```

## Steady state solution
For steady state flow solutions, equations {eq}`eqn:eta_HalfPlaneModel` to {eq}`eqn:p_HalfPlaneModel` can be further simplified since $\omega=0$ and hence $\Omega=-u_0k$. For completeness, the simplified solution is given below:
```{math}
\hat{\eta}_1(k,z) = \hat{h}(k)e^{jm(k)z},
```
with
```{math}
:label: eqn:h_k
    \hat{h}(k)=\frac{1}{2\pi}\int_{-\infty}^{\infty}h(x)\exp{(-jkx)}\;\mathrm{d} x
```
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