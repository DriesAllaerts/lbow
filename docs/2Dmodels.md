(sec:2DModels)=
# Two-dimensional models

Two-dimensional models allow the flow perturbation to vary in two direction, like an isolated mountain.


(sec:PolarizationEqns)=
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