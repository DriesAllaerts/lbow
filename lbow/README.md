# Linear Buoyancy Wave package (LBoW)
LBoW is a python package to solve Linear Buoyancy Wave problems. The following modules and classes are available:
- [steadystate](steadystate.py) module, containing models to compute the steady state flow solution
  - `OneLayerModel` class: Base class for steady state models consisting of one single layer (i.e., uniform background flow parameters)
    - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
    See [steady-state--1-surface-corrugation](../notebooks/steady-state--1-surface-corrugation.ipynb) and [steady-state--2-witch-of-Agnesi](../notebooks/steady-state--2-witch-of-Agnesi.ipynb) notebooks for usage.
    - `ChannelModel` class: Model with the ground surface as bottom boundary and a rigid lid as the top boundary.  
    See [steady-state--3-trapped-wave-solution](../notebooks/steady-state--3-trapped-wave-solution.ipynb) notebook for usage.
- [transient](transient.py) module, containing models to compute transient flow behaviour
  - `OneLayerModel` class: Base class for transient models consisting of one single layer (i.e., uniform background flow parameters)
  - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition).  
  See [transient--1-impulse-response](../notebooks/transient--1-impulse-response.ipynb) and [transient--2-step-function](../notebooks/transient--2-step-function.ipynb) notebooks for usage.
  
Currently, only models for one-dimensional perturbations (e.g., ridge lines) leading to two-dimensional flow problems (x-z) are implemented.
