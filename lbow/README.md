# Linear Buoyancy Wave package (LBoW)
LBoW is a python package to solve Linear Buoyancy Wave problems. The following modules and classes are available:
- `steadystate` module, containing models to compute the steady state flow solution
  - `OneLayerModel` class: Base class for steady state models consisting of one single layer (i.e., uniform background flow parameters)
  - `ChannelModel` class: Model with the ground surface as bottom boundary and a rigid lid as the top boundary
  - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition)
- `transient` module, containing models to compute transient flow behaviour
  - `OneLayerModel` class: Base class for transient models consisting of one single layer (i.e., uniform background flow parameters)
  - `HalfPlaneModel` class: Model with the ground surface as bottom boundary and the top boundary at infinity (radiation boundary condition)
  
Currently, only models for one-dimensional perturbations (e.g., ridge lines) leading to two-dimensional flow problems (x-z) are implemented.
