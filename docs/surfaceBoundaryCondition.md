(sec:surfaceBC)=
# Surface boundary condition
Topographical features are usually included in linear theory by requiring that the vertical displacement of the flow at the surface matches the surface elevation {cite}`smith_linear_1980`, giving rise to the linearized boundary condition

$$
\eta_1(x,y,0,t) = h(x,y,t).
$$ (eqn:real_BC)

```{note}
For the time scales relevant to wind flow over hills and mountains, the shape of the topography is normally not time dependent. However, time is kept here as an independent variable to allow the study of theoretical cases where the surface perturbation is time-dependent. Such a case is for example useful to visualize how buoyancy waves start to develop in numerical simulations before a steady state is reached.
```
In general, $h(x,y,t)$ will excite a spectrum of plane waves, so the solution to equation {eq}`eqn:TaylorGoldstein_eta` will be of the form

$$
\eta_1(x,y,z,t)=\iiint_{-\infty}^{\infty}\hat{\eta}_1(k,l,z,\omega)\exp{[j(kx+ly-\omega t)]}\;\mathrm{d} k \, \mathrm{d} l \, \mathrm{d} \omega
$$

To obtain a boundary condition for $\hat{\eta}_1(k,l,\omega,z)$, the surface elevation $h(x,y,t)$ is decomposed into planar waves as well (in other words, the Fourier transform is applied)

$$

    \hat{h}(k,l,\omega)=\frac{1}{8\pi^3}\iiint_{-\infty}^{\infty}h(x,y,t)\exp{[-j(kx+ly-\omega t)]}\;\mathrm{d} x \,\mathrm{d} y \,\mathrm{d} t
$$

The boundary condition in real space (equation {eq}`eqn:real_BC`) then translates into a boundary condition per wavenumber--frequency pair:

$$
    \hat{\eta}_1(k,l,0,\omega) = \hat{h}(k,l,\omega).
$$ (eqn:Fourier_BC)

%Next to the surface boundary condition, an upper boundary condition is needed. Two cases are considered in sections {numref}`%s <sec:halfPlaneModel>` and {numref}`%s <sec:ChannelModel>`