(sec:2DModels)=
# Two-dimensional models

Two-dimensional models allow the flow perturbation to vary in two direction, like an isolated mountain.


(sec:PolarizationEqns2D)=
## Polarization equations

As with one-dimensional models, flow variables like $u_1$, $v_1$, $w_1$, and $p_1$ are easily related to the vertical displacement $\hat{\eta}_1(k,l,z,\omega)$. As before, the vertical velocity $w_1$ follows from the kinematic condition {eq}`eqn:kinematicCondition`:

$$
    \hat{w}_1(k,l,z,\omega)=-j\Omega\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:polarization_w_2D)

We can no longer obtain the streamwise velocity component $u_1$ from the continuity equation {eq}`eqn:cont` as in 2D the continuity equation also depends on the spanwise velocity component $v_1$. Instead, we first need to find the pressure $p_1$. The presure follows from equation {eq}`eqn:horizontal_part`

$$
    \hat{p}_1(k,l,z,\omega)=\rho_0\frac{j\Omega}{k^2+l^2}\frac{\partial}{\partial z}\hat{w}_1(k,l,z,\omega) = \rho_0\frac{\Omega^2}{k^2+l^2}\frac{\partial}{\partial z}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:polarization_p_2D)

The horizontal velocity components then follow from the momentum equations {eq}`eqn:umom`-{eq}`eqn:vmom`:

$$
    \hat{u}_1(k,l,z,\omega) = \frac{1}{\rho_0}\frac{k}{\Omega} \hat{p}_1(k,l,z,\omega) = \frac{\Omega k}{k^2+l^2}\frac{\partial}{\partial z}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:polarization_u_2D)

$$
    \hat{v}_1(k,l,z,\omega) = \frac{1}{\rho_0}\frac{l}{\Omega} \hat{p}_1(k,l,z,\omega) = \frac{\Omega l}{k^2+l^2}\frac{\partial}{\partial z}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:polarization_v_2D)

(sec:HalfPlaneModel2D)=
## Half Plane Model

As in the 1D case, the 2D Half Plane Model solves the flow on an infinite half-plane above the ground. The structure of this section follows that of {numref}`Section %s <sec:HalfPlaneModel>`.

### Transient solution

Waves are only excited at the surface, and therefore the solution only contains waves that propagate/radiate upwards. Hence, we retain the term with $m_1=m$ in the general solution {eq}`eqn:generalSolution` and discard the term with $m_2=-m$. This gives

$$
\hat{\eta}_1(k,l,z,\omega) = A(k,l,\omega)e^{jm(k,l,\omega)z},
$$ (eqn:eta_HalfPlaneModel2D)

with $A(k,l,\omega)=\hat{h}(k,l,\omega)$ following from equation {eq}`eqn:Fourier_BC`, and the vertical wave number defined in equation {eq}`eqn:verticalWaveNumber`.

With the vertical displacement given by equations {eq}`eqn:eta_HalfPlaneModel2D` and {eq}`eqn:verticalWaveNumber`, the other flow variables can be calculated using the {ref}`polarization equations <sec:polarizationEqns2D>`:

$$
    \hat{w}_1(k,l,z,\omega)=-j\Omega\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:w_HalfPlaneModel2D)

$$
    \hat{p}_1(k,l,z,\omega)=\rho_0\frac{j\Omega^2 m}{k^2+l^2}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:p_HalfPlaneModel2D)

$$
    \hat{u}_1(k,l,z,\omega)=\frac{j\Omega km}{k^2+l^2}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:u_HalfPlaneModel2D)

$$
    \hat{v}_1(k,l,z,\omega)=\frac{j\Omega lm}{k^2+l^2}\hat{\eta}_1(k,l,z,\omega)
$$ (eqn:v_HalfPlaneModel2D)



The full flow solution for any flow variable is found by applying the 3D inverse Fourier transform:

$$
    \phi_1(x,y,z,t)=\iiint_{-\infty}^{\infty}\hat{\phi}_1(k,l,z,\omega)\exp{[j(kx+ly-\omega t)]}\;\mathrm{d} k \,\mathrm{d} l \,\mathrm{d} \omega
$$

with $\phi_1$ representing any variable in $(\eta_1,u_1,v_1,w_1,p_1)$.


### Steady state solution
For steady state flow solutions, $\omega=0$ and hence $\Omega=-u_0k-v_0l$. Given this form of the intrinsic frequency, the flow solution is given by the same set of equations as in the 2D transient case except now the solution no longer depends on $\omega$ (and hence the inverse Fourier transform is only over two instead of three dimensions).