# Linear Buoyancy Wave package (LBoW)
LBoW is a python package to solve Linear Buoyancy Wave problems. The following modules and classes are available:
- oneDimensional module, containing models for one-dimensional perturbations (ridge lines) leading to two-dimensional flow problems (x-z)
  - [steadystate](steadystate.py) submodule, containing models to compute the steady state flow solution
    - `OneLayerModel` class: Base class for steady state models consisting of one single layer (i.e., uniform background flow parameters)
    - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
    See [1-1D-steady-state--surface-corrugation](../notebooks/1-1D-steady-state--surface-corrugation.ipynb) and [2-1D-steady-state--witch-of-Agnesi](../notebooks/2-1D-steady-state--witch-of-Agnesi.ipynb) notebooks for usage.
    - `ChannelModel` class: Model with the ground surface as bottom boundary and a rigid lid as the top boundary.  
    See [3-1D-steady-state--trapped-wave-solution](../notebooks/3-1D-steady-state--trapped-wave-solution.ipynb) notebook for usage.
  - [transient](transient.py) submodule, containing models to compute transient flow behaviour
    - `OneLayerModel` class: Base class for transient models consisting of one single layer (i.e., uniform background flow parameters)
    - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
    See [4-1D-transient--impulse-response](../notebooks/4-1D-transient--impulse-response.ipynb) and [5-1D-transient--step-function](../notebooks/5-1D-transient--step-function.ipynb) notebooks for usage.
- twoDimensional module, containing models for two-dimensional perturbations (isolated mountains) leading to three-dimensional flow problems (x-y-z)
