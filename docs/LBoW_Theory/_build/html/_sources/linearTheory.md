(sec:LinearTheory)=
# Linear theory -- General form

It is assumed that the velocity $(u_0,v_0,0)$ and the Brunt--Väisälä or buoyancy frequency $N$ of the background state are independent of height, and that the Boussinesq approximation can be made. In the absence of rotation and friction, the governing equations for small perturbations $(u_1,v_1,w_1)$ in velocity, $p_1$ in pressure and $\theta_1$ in potential temperature are
```{math}
:label: eqn:umom
   \frac{\mathrm{D}u_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial x},
```
```{math}
:label: eqn:vmom
   \frac{\mathrm{D}v_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial y},
```
```{math}
:label: eqn:wmom
    \frac{\mathrm{D}w_1}{\mathrm{D}t} =-\frac{1}{\rho_0}\frac{\partial p_1}{\partial z} + \frac{\theta_1}{\theta_0}g,
```
```{math}
:label: eqn:theta
    \frac{\mathrm{D}\theta_1}{\mathrm{D}t} + w_1\frac{\mathrm{d}\theta_0}{\mathrm{d}z} =0,
```
```{math}
:label: eqn:cont
    \frac{\partial u_1}{\partial x}+\frac{\partial v_1}{\partial y}+\frac{\partial w_1}{\partial z}=0
```
with the material derivative $\frac{\mathrm{D}}{\mathrm{D}t} = \left(\frac{\partial}{\partial t}+\boldsymbol{u}_0\cdot\boldsymbol{\nabla}\right)$. Equations {eq}`eqn:umom`-{eq}`eqn:cont` can be reduced to a single equation in $w_1$ as follows. Taking the material derivative of the vertical momentum equation {eq}`eqn:wmom` allows the substitution of the potential temperature equation {eq}`eqn:theta`. The pressure is found by taking the divergence of the momentum equations {eq}`eqn:umom`-{eq}`eqn:wmom` and applying the continuity equation {eq}`eqn:cont`. This yields
```{math}
:label: eqn:TaylorGoldstein_w
    \left(\frac{\mathrm{D}}{\mathrm{D}t}\right)^2\nabla^2 w_1+N^2\nabla^2_H w_1=0.
```
with the horizontal Laplacian operator and Brunt--Väisälä frequency defined as $\nabla^2_H=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}$ and $N^2=\frac{g}{\theta_0}\frac{\mathrm{d}\theta_0}{\mathrm{d}z}$, respectively. Equation~{eq}`eqn:TaylorGoldstein_w` is a simplified form of the more general Taylor-Goldstein equation for wave motions in a stably stratified shear flow.

Equation {eq}`eqn:TaylorGoldstein_w` can also be expressed in terms of the vertical displacement $\eta_1$ of a fluid parcel above its undisturbed level, which is related to the vertical wind speed perturbation via the kinematic condition {cite}`smith_linear_1980`
```{math}
:label: eqn:kinematicCondition
    w_1=\frac{\mathrm{D}\eta_1}{\mathrm{D}t}.
```
Substituting equation {eq}`eqn:kinematicCondition` into equation {eq}`eqn:TaylorGoldstein_w` gives
```{math}
:label: eqn:TaylorGoldstein_eta
    \left(\frac{\mathrm{D}}{\mathrm{D}t}\right)^2\nabla^2 \eta_1+N^2\nabla^2_H \eta_1=0.
```
Assuming a plane wave solution
```{math}
    \eta_1(x,y,z,t)=\hat{\eta}_1(z)\exp{[j(kx+ly-\omega t)]}
```
results in
```{math}
:label: eqn:characteristic_eqn
    \frac{\partial^2\hat{\eta}_1}{\partial z^2}+m^2\hat{\eta}_1=0
```
with the vertical wave number $m$ given by
results in
```{math}
:label: eqn:dispersionEquation
    m^2=(k^2+l^2)\left(\frac{N^2}{\Omega^2}-1\right).
```
The intrinsic frequency is hereby defined as $\Omega=\omega-\mathbf{u_0}\cdot\mathbf{k}$, with $\mathbf{k}=(k,l,m)$ the wave vector and $\mathbf{u_0}=(u_0,v_0,0)$ the background wind speed vector. It is important to remark that equation {eq}`eqn:dispersionEquation` gives rise to two possible solutions for $m$. For $\Omega^2<N^2$, $m$ is a real number and the planar wave is vertically propagating, while for $\Omega^2>N^2$, $m$ is imaginary and the wave becomes evanescent.

Equation {eq}`eqn:characteristic_eqn` is a linear, homogeneous, ordinary differential equation of second order, for which the general solution can be written as
```{math}
\hat{\eta}_1(k,l,z,\omega) = A(k,l,\omega)e^{jm_1(k,l,\omega)z}+B(k,l,\omega)e^{jm_2(k,l,\omega)z},
```
with $m_1$ and $m_2$ the positive and negative roots of equation {eq}`eqn:dispersionEquation`, respectively. The coefficients $A$ and $B$ are determined by the boundary conditions.