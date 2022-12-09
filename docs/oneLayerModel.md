(sec:OneLayerModel)=
# One-Layer Models (1D)

In this section, the solution to equation {eq}`eqn:characteristic_eqn` is derived for various one layer models, i.e., models consisting of one layer with uniform background wind speed and buoyancy frequency. The flow perturbation is triggered by the bottom boundary condition, which represents a topographical feature. The difference between the models stems from the upper boundary condition. All models assume the flow perturbation is one-dimensional, like a ridge line, for which there is no variation in the spanwise direction and $l=0$.

Topographical features are usually included in linear theory by requiring that the vertical displacement of the flow at the surface matches the surface elevation {cite}`smith_linear_1980`, giving rise to the linearized boundary condition for one-dimensional topographical features
```{math}
:label: eqn:real_BC
\eta_1(x,0,t) = h(x,t).
```
```{note}
For the time scales relevant to wind flow over hills and mountains, the shape of the topography is normally not time dependent. However, time is kept here as an independent variable to allow the study of theoretical cases where the surface perturbation is time-dependent. Such a case is for example useful to visualize how buoyancy waves start to develop in numerical simulations before a steady state is reached.
```
In general, $h(x,t)$ will excite a spectrum of plane waves, so the solution to equation {eq}`eqn:TaylorGoldstein_eta` will be of the form
```{math}
\eta_1(x,z,t)=\iint_{-\infty}^{\infty}\hat{\eta}_1(k,z,\omega)\exp{[j(kx-\omega t)]}\;\mathrm{d} k \,\mathrm{d} \omega
```
To obtain a boundary condition for $\hat{\eta}_1(k,\omega,z)$, the surface elevation $h(x,t)$ is decomposed into planar waves as well (in other words, the two-dimensional Fourier transform is applied)
```{math}
    \hat{h}(k,\omega)=\frac{1}{4\pi^2}\iint_{-\infty}^{\infty}h(x,t)\exp{[-j(kx-\omega t)]}\;\mathrm{d} x \,\mathrm{d} t
```
The boundary condition in real space (equation {eq}`eqn:real_BC`) then translates into a boundary condition per wavenumber--frequency pair:
```{math}
:label: eqn:Fourier_BC
    \hat{\eta}_1(k,0,\omega) = \hat{h}(k,\omega).
```

Next to the surface boundary condition, an upper boundary condition is needed. Two cases are considered in the following subsections.