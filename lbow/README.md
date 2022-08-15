# Linear Buoyancy Wave package (LBoW)
LBoW is a python package to solve Linear Buoyancy Wave problems. The following modules and classes are available:
- `steadystate` module, containing models to compute the steady state flow solution
  - `OneLayerModel` class: Base class for steady state models consisting of one single layer (i.e., uniform background flow parameters)
  - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
  See [1-steady-state--surface-corrugation](../notebooks/1-steady-state--surface-corrugation.ipynb) and [2-steady-state--witch-of-Agnesi](../notebooks/2-steady-state--witch-of-Agnesi.ipynb) notebooks for usage.
  - `ChannelModel` class: Model with the ground surface as bottom boundary and a rigid lid as the top boundary.  
  See [3-steady-state--trapped-wave-solution](../notebooks/3-steady-state--trapped-wave-solution.ipynb) notebook for usage.
- `transient` module, containing models to compute transient flow behaviour
  - `OneLayerModel` class: Base class for transient models consisting of one single layer (i.e., uniform background flow parameters)
  - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
  See [4-transient--impulse-response](../notebooks/4-transient--impulse-response.ipynb) and [5-transient--step-function](../notebooks/5-transient--step-function.ipynb) notebooks for usage.
  
Currently, only models for one-dimensional perturbations (e.g., ridge lines) leading to two-dimensional flow problems (x-z) are implemented.
