(sec:MultiLayerModel)=
# Multi-Layer Models (1D)

In this section, the solution to equation {eq}`eqn:characteristic_eqn` is derived for a multi-layer model, i.e., a model consisting of multiple layers each with uniform background wind speed and buoyancy frequency. The flow perturbation is triggered by the bottom boundary condition in the bottom layer adjacent to the surface. The upper layer is assumed to be semi-infinite, and the top boundary conditions for the upper layer is a radiation condition. At the interface between different layers, the flow solution should be continuous, and this requirement leads to interface conditions. All models assume the flow perturbation is one-dimensional, like a ridge line, for which there is no variation in the spanwise direction and $l=0$.

For now, we only consider a steady state model.

## Steady state solution

Assume we want to build a steady state model with $Q$ layers ($q=0,1,\ldots,Q-1$). The layers are defined by an array of heights $z_q$ such that layer $q$ extents from $z_q$ to $z_{q+1}$, layer $q+1$ extents from $z_{q+1}$ to $z_{q+2}$, etc. The velocity and buoyancy frequency in layer $q$ are denoted as $U_q$ and $N_q$. Within each layer, we need to solve equation {eq}`eqn:characteristic_eqn` for the displacement $\eta_q$, with the vertical wave number $m_q$ computed with equation {eq}`eqn:dispersionEquation` based on $U_q$ and $N_q$. The general solution in each layer can be written as

$$
\hat{\eta}_q(k,z)=A_q e^{jm_q(z-z_q)}+B_q e^{jm_q(z_{q+1}-z)}.
$$ (eqn:characteristic_eqn_layerq)

The sign of the vertical wave number $m_q$ is chosen such that $A_q e^{jm_q(z-z_q)}$ is propagating upward and $B_q e^{jm_q(z_{q+1}-z)}$ is propagating downward (same sign convention as for the {ref}`one-layer half-plane model <sec:HalfPlaneModel>`).

```{note}
We use the base functions $A_q e^{jm_q(z-z_q)}$ and $B_q e^{jm_q(z_{q+1}-z)}$ rather than $A_q e^{jm_qz}$ and $B_q e^{-jm_qz}$ to deal with overflow in the exponential for evanescent waves. For evanescent waves, propagating upward means the amplitude decreases as the height increases, whereas propagating downward means the amplitude decreases as the height decreases. The height offset for upward propagating waves is chosen such that the amplitude is equal to $A_q$ at the offset $z_q$ (corresponding to the base height of layer $q$), whereas for downward propagating waves the amplitude is equal to $B_q$ at $z_{q+1}$ (corresponding to the top height of layer $q$). In other words, within layer $q$, the maximum amplitude of upward travelling evanescent waves is $A_q$ and occurs at the bottom of the layer ($z_q$), while the maximum aplitude of downward travelling evanescent waves is $B_q$ and occurs at the top of the layer ($z_{q+1}$).
```

The surface boundary conditions $\eta(x,0)=h(x)$ implies

$$
A_0+B_0 e^{jm_0(z_1-z_0)} = \hat{h}.
$$ (eqn:surface_BC_multiLayer)

The requirement that the solution remains bounded for $z\rightarrow+\infty$ implies

$$
B_{Q-1}=0
$$ (eqn:top_BC_multiLayer)

The requirement of the solution being continuous at the interfaces can be expressed as

$$
\eta_{q-1}(z_q)=\eta_q(z_q)
$$

$$
\frac{\partial \eta_{q-1}}{\partial z}(z_q)=\frac{\partial \eta_q}{\partial z}(z_q)
$$

for $q=1,\ldots,Q-1$.
```{note}
In some studies {cite}`vosper_inversion_2004,devesse_including_2022` interface conditions are expressed in terms of continuity of the vertical velocity and the pressure. I'm not sure whether this is the same as saying that $\eta$ and $\partial_z \eta$ should be continuous. To be further investigated.
```

Filling in equation {eq}`eqn:characteristic_eqn_layerq` in the interface conditions leads to

$$
A_{q-1}e^{jm_{q-1}\Delta z_{q-1}}+B_{q-1}=A_q+B_qe^{jm_q\Delta z_q)}
$$ (eqn:interface_eta_BC)

$$
m_{q-1}\left(A_{q-1}e^{jm_{q-1}\Delta z_{q-1}}-B_{q-1}\right)=m_q\left(A_q-B_qe^{jm_q\Delta z_q}\right)
$$ (eqn:interface_deta_BC)

for $q=1,\ldots,Q-1$. The thickness of layer $q$ is denoted as $\Delta z_q=z_{q+1}-z_q$. Equations {eq}`eqn:surface_BC_multiLayer`-{eq}`eqn:interface_deta_BC` combined result in a system of $2Q$ equations in $2Q$ unknown coefficients $A_q$ and $B_q$. This system of equations depends on $k$ and needs to be solved for every value of $k$.
```{note}
For the zero mode $k=0$, equation {eq}`eqn:dispersionEquation` leads to $m=0$ and the system of boundary conditions becomes singular. The boundary conditions only say that $A_0+B_0=A_1+B_1=\ldots=\hat{h}$. Hence, for $k=0$, we choose to set $A_0=A_1=\ldots=\hat{h}$ and $B_0=B_1=\ldots=0$.
```

Once $A_q$ and $B_q$ have been found, the displacement (in Fourier space) is given by equation {eq}`eqn:characteristic_eqn_layerq`. The vertical velocity is given by 

$$
\hat{w}_q(k,z)=jU_qkA_q e^{jm_q(z-z_q)}+jU_qkB_q e^{jm_q(z_{q+1}-z)}.
$$

The full solution for a perturbation quantity $\hat{\phi}$ is then given by a piecewise function

$$
    \hat{\phi}(k,z) = \begin{cases}
    \hat{\phi}_0(k,z) & \text{for}\;z_0\le z<z_1, \\
    \hat{\phi}_1(k,z) & \text{for}\;z_1\le z<z_2, \\
    \ldots \\
    \hat{\phi}_{Q-1}(k,z) & \text{for}\;z_{Q-1}\le z
    \end{cases}
$$

with $\hat{\phi}$ representing any variable in $(\hat{\eta},\hat{u},\hat{w},\hat{p})$
