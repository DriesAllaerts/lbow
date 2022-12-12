(sec:ChannelModel)=
# Channel Model

The Channel Model assumes that the flow is bounded above by a rigid lid at a certain height $H$, causing downward wave reflections. This situation is reminiscent of numerical simulations of buoyancy waves with inadequate reflective boundary conditions at the top of the numerical domain.

We first derive the solution for the general, transient flow scenario. Next, we present the simpler solution for a steady state flow scenario.

## Transient solution

The rigid lid condition translates to $\eta_1(x,H,t)=0$, and hence $\hat{\eta}_1(k,H,\omega)=0$. Combined with the surface boundary condition {eq}`eqn:Fourier_BC`, this results in a system of equations in terms of the unknowns $A$ and $B$:

$$
	\begin{cases}
		A+B &= \hat{h}(k,\omega) \\
		Ae^{jmH}+Be^{-jmH} &= 0
	\end{cases}
$$ (eqn:ABsystem)

where we made use of the sign convention for the vertical wave number defined in {numref}`Section %s <sec:verticalWaveNumber>`. The solution to {eq}`eqn:ABsystem` can be written as

$$
	A=\hat{h}(k,\omega)/(1-e^{2jmH}) \quad\text{and}\quad B = -\hat{h}(k,\omega)e^{2jmH}/(1-e^{2jmH})
$$

The flow solution is then given by

$$
	\hat{\eta}_1(k,z,\omega)=\frac{e^{jmz}-e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k,\omega),
$$

and

$$
    \hat{w}_1(k,z,\omega)=-j\Omega\hat{\eta}_1(k,z,\omega)=-j\Omega\frac{e^{jmz}-e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k,\omega),
$$

$$
	\hat{u}_1(k,z,\omega)=-\frac{1}{jk}\frac{\partial}{\partial z}\hat{w}_1(k,z,\omega)=\frac{jm\Omega}{k}\frac{e^{jmz}+e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k,\omega),
$$

$$
	\hat{p}_1(k,z,\omega)=\rho_0 \frac{\Omega}{k}\hat{u}_1(k,z,\omega)=\rho_0\frac{jm\Omega^2}{k^2}\frac{e^{jmz}+e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k,\omega).
$$

## Steady state solution
As before, the general solution can be simplified under steady state conditions since $\omega=0$ and hence $\Omega=-u_0k$. The solution is

$$
	\hat{\eta}_1(k,z)=\frac{e^{jmz}-e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k),
$$

with $\hat{h}(k)$ and $m(k)$ defined by equations {eq}`eqn:h_k` and {eq}`eqn:m_k`, respectively. The solution for other flow variables is given by

$$
    \hat{w}_1(k,z)=ju_0k\hat{\eta}_1(k,z)=ju_0k\frac{e^{jmz}-e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k),
$$

$$
	\hat{u}_1(k,z)=-\frac{1}{jk}\frac{\partial}{\partial z}\hat{w}_1(k,z)=-jmu_0\frac{e^{jmz}+e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k),
$$

$$
	\hat{p}_1(k,z)=-\rho_0 u_0\hat{u}_1(k,z)=\rho_0 jmu_0^2\frac{e^{jmz}+e^{jm(2H-z)}}{(1-e^{2jmH})}\hat{h}(k).
$$