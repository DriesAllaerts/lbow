(sec:verticalWaveNumber)=
# Vertical wave number

As shown previously, the vertical wave number of bouyancy waves is related to the buoyancy frequency, the intrinsic frequency, and the horizontal wave numbers by means of the dispersion equation:

$$
    m^2=(k^2+l^2)\left(\frac{N^2}{\Omega^2}-1\right).
$$ (eqn:dispersionEquation)

This equation gives rise to two possible wave types. For $\Omega^2<N^2$, $m$ is a real number and the planar wave will be vertically propagating, while for $\Omega^2>N^2$, $m$ is imaginary and the wave becomes evanescent.

For both propagating and evanescent waves, there are two roots to equation {eq}`eqn:dispersionEquation` and hence two terms in the general solution {eq}`eqn:generalSolution`. For propagating waves, the two terms correspond to upward and downward propagating waves. For evanescent waves,  the two terms in the general solution differ in terms of the direction in which the wave amplitude attenuates. One term corresponds to waves that decrease in amplitude as the height *increases*, while the other term corresponds to waves that decrease in amplitude as the height *decreases*. An example of the former is an evanescent wave triggered at the surface that radiates upwards, while the latter could be an evanescent wave radiating downwards due to a reflection of an upward radiating wave at a certain height.

It is convenient for the further development of the theory to introduce a sign convention for the vertical wave number $m$, so that the roots to equation {eq}`eqn:dispersionEquation` are given by $m_1=m$ and $m_2=-m$. We choose $m$ to correspond to propagating/evanescent waves that propagate/radiate upwards. Downward propagating/radiating waves are then represented by $-m$.

For evanescent waves ($\Omega^2>N^2$), the general solution {eq}`eqn:generalSolution` shows that we need to take the positive imaginary root to obtain a decreasing amplitude with increasing height.

For propagating waves ($\Omega^2<N^2$), we need to identify waves with a vertical group velocity pointing upwards. The group velocity relative to the flow is defined as $\mathbf{c}_g=(\partial\Omega/\partial k,\partial\Omega/\partial l,\partial\Omega/\partial m)$. The group velocity relative to the ground can be found by taking the derivative of the apparent frequency $\omega$ rather than the intrinsic frequency, but since $\omega=\Omega+\mathbf{u}_0\cdot\mathbf{k}$, this comes down to a simple vector summation. As the background vertical velocity is assumed to be zero, $w_g=\partial\Omega/\partial m$. The partial derivative of $\Omega$ to $m$ can be found by reformulating the expression of the vertical wavenumber $m$ ({eq}`eqn:dispersionEquation`) in terms of $\Omega$, which gives $\Omega=\pm N\sqrt{k^2+l^2}/\sqrt{k^2+l^2+m^2}$. Then, it can be shown that $w_g=-\Omega m/(k^2+l^2+m^2)$. Hence, for $w_g$ to be positive, the vertical wave number must be chosen so that $\text{sign}(m)=-\text{sign}(\Omega)$.

In summary, upward propagating/radiating waves have a vertical wave number given by

$$
    m(k,l,\omega) = \begin{cases}
    j\sqrt{k^2+l^2}\sqrt{1-N^2/\Omega^2} & \text{for}\;\Omega^2>N^2, \\
    -\text{sign}(\Omega)\,\sqrt{k^2+l^2}\sqrt{N^2/\Omega^2-1} & \text{for}\;\Omega^2<N^2.
    \end{cases}
$$ (eqn:verticalWaveNumber)

For downward propagating/radiating waves, the vertical wave number is $-m$.

(sec:hydrostaticAssumtion)=
## Hydrostatic assumption
Sometimes the flow perturbation is assumed to be in hydrostatic equilibrium. Mathematically, this means that in equation {eq}`eqn:wmom`, the left-hand side is assumed to be negligible and the vertical pressure gradient is balanced by the buoyancy force. In this case, governing equations can still be reduced to equation {eq}`eqn:characteristic_eqn` but with the following dispersion relation:

$$
    m^2=(k^2+l^2)\frac{N^2}{\Omega^2}.
$$ (eqn:hydrostaticDispersionEquation)

The hydrostatic dispersion relation will always yield a real number for the vertical wave number, which means that under the hydrostatic assumption all planar waves propagate vertically and there are no evanescent waves.

The vertical group velocity can be found in the same way as before. The intrinsic frequency is now given by $\Omega=\pm N\sqrt{k^2+l^2}/|m|$, from which it follows that $w_g=\partial\Omega/\partial m=-\Omega /m$. Similar to the non-hydrostatic case, for $w_g$ to be positive, the vertical wave number must be chosen so that $\text{sign}(m)=-\text{sign}(\Omega)$.

In summary, upward propagating hydrostatic waves have a vertical wave number given by

$$
    m(k,l,\omega) = 
    -\sqrt{k^2+l^2}\,\frac{N}{\Omega} 
$$ (eqn:verticalWaveNumber)