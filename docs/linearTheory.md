(sec:LinearTheory)=
# Linear theory -- General form

It is assumed that the velocity $(u_0,v_0,0)$ and the Brunt--Väisälä or buoyancy frequency $N$ of the background state are independent of height, and that the Boussinesq approximation can be made. In the absence of rotation and friction, the governing equations for small perturbations $(u_1,v_1,w_1)$ in velocity, $p_1$ in pressure and $\theta_1$ in potential temperature are

$$
    \frac{\mathrm{D}u_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial x},
$$ (eqn:umom)

$$
   \frac{\mathrm{D}v_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial y},
$$ (eqn:vmom)

$$
    \frac{\mathrm{D}w_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial z} + \frac{\theta_1}{\theta_0}g,
$$ (eqn:wmom)

$$
    \frac{\mathrm{D}\theta_1}{\mathrm{D}t} + w_1\frac{\mathrm{d}\theta_0}{\mathrm{d}z} =0,
$$ (eqn:theta)

$$
    \frac{\partial u_1}{\partial x}+\frac{\partial v_1}{\partial y}+\frac{\partial w_1}{\partial z}=0
$$ (eqn:cont)

with the material derivative $\frac{\mathrm{D}}{\mathrm{D}t} = \left(\frac{\partial}{\partial t}+\boldsymbol{u}_0\cdot\boldsymbol{\nabla}\right)$. Equations {eq}`eqn:umom`-{eq}`eqn:cont` can be reduced to a single equation in $w_1$ as follows {cite}`gill_atmosphere-ocean_1982`. Taking the material derivative of the vertical momentum equation {eq}`eqn:wmom` allows the substitution of the potential temperature equation {eq}`eqn:theta`:

$$
    \frac{\mathrm{D}^2 w_1}{\mathrm{D}t^2} +N^2 w_1 = -\frac{1}{\rho_0}\frac{\mathrm{D}}{\mathrm{D}t}\frac{\partial p_1}{\partial z}
$$ (eqn:vertical_part)

with the buoyancy frequency defined as $N^2=\frac{g}{\theta_0}\frac{\mathrm{d}\theta_0}{\mathrm{d}z}$. An equation for the pressure can be found by taking the time derivative of the continuity equation {eq}`eqn:cont` and substituting the momentum equations {eq}`eqn:umom`-{eq}`eqn:vmom`:

$$
    \frac{\mathrm{D}}{\mathrm{D}t} \frac{\partial w_1}{\partial z} = \frac{1}{\rho_0}\nabla^2_H p_1.
$$ (eqn:horizontal_part)

with the horizontal Laplacian operator defined as $\nabla^2_H=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}$. Combining equations {eq}`eqn:vertical_part` and {eq}`eqn:horizontal_part` to eliminate $p_1$ yields

$$
    \left(\frac{\mathrm{D}}{\mathrm{D}t}\right)^2\nabla^2 w_1+N^2\nabla^2_H w_1=0.
$$ (eqn:TaylorGoldstein_w)

Equation {eq}`eqn:TaylorGoldstein_w` is a simplified form of the more general Taylor-Goldstein equation for wave motions in a stably stratified shear flow.

Equation {eq}`eqn:TaylorGoldstein_w` can also be expressed in terms of the vertical displacement $\eta_1$ of a fluid parcel above its undisturbed level, which is related to the vertical wind speed perturbation via the kinematic condition {cite}`smith_linear_1980`

$$
    w_1=\frac{\mathrm{D}\eta_1}{\mathrm{D}t}.
$$ (eqn:kinematicCondition)

Substituting equation {eq}`eqn:kinematicCondition` into equation {eq}`eqn:TaylorGoldstein_w` gives

$$
    \left(\frac{\mathrm{D}}{\mathrm{D}t}\right)^2\nabla^2 \eta_1+N^2\nabla^2_H \eta_1=0.
$$ (eqn:TaylorGoldstein_eta)

We can find a solution for equation {eq}`eqn:TaylorGoldstein_eta` by assuming a plane wave solution of the form

$$
    \eta_1(x,y,z,t)=\hat{\eta}_1(z)\exp{[j(kx+ly-\omega t)]}.
$$

Substituting the planar wave form into {eq}`eqn:TaylorGoldstein_eta`
results in

$$
    \frac{\partial^2\hat{\eta}_1}{\partial z^2}+m^2\hat{\eta}_1=0
$$ (eqn:characteristic_eqn)

with the vertical wave number $m$ given by

$$
    m^2=(k^2+l^2)\left(\frac{N^2}{\Omega^2}-1\right).
$$ (eqn:dispersionEquation_1)

The intrinsic frequency is hereby defined as $\Omega=\omega-\mathbf{u_0}\cdot\mathbf{k}$, with $\mathbf{k}=(k,l,m)$ the wave vector and $\mathbf{u_0}=(u_0,v_0,0)$ the background wind speed vector.

Equation {eq}`eqn:characteristic_eqn` is a linear, homogeneous, ordinary differential equation of second order, for which the general solution can be written as

$$
\hat{\eta}_1(k,l,z,\omega) = A(k,l,\omega)e^{jm_1(k,l,\omega)z}+B(k,l,\omega)e^{jm_2(k,l,\omega)z},
$$ (eqn:generalSolution)

with $m_1$ and $m_2$ the roots of equation {eq}`eqn:dispersionEquation_1`. The coefficients $A$ and $B$ are determined by the boundary conditions.