# Linear Buoyancy Wave package (LBoW)
LBoW is a python package to solve Linear Buoyancy Wave problems.

See the `dev` branch for the latest code under development

## Requirements
The code runs from python 3.10, but highest version is recommended. The core python libraries of the package only depend on numpy. Demonstration of the package is done via jupyter notebooks and depends on jupyter and matplotlib.

## Installation
To install, run `pip install -e lbow` after cloning the repository (or `pip install -e .` from inside `/path/to/lbow`).

## Structure
In this repository you can find:
- `lbow/` directory containing the core python libraries
- `notebooks/` directory containing jupyter notebooks to demonstrate the use of the lbow package
- `docs/` directory containing documentation regarding the theory behind this package (see also below)

## Documentation
The theory and mathematical expressions behind LBoW are described in a Jupyter Book that is hosted as a GitHub Page at [https://driesallaerts.github.io/lbow/intro.html](https://driesallaerts.github.io/lbow/intro.html).

The source code for the Jupyter Book can be found in the `docs/` directory. The book's build files (i.e., the static html files) live in the `gh-pages` branch.

## Author(s)
This software has been developed by
**Dries Allaerts** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0002-8758-3952](https://orcid.org/0000-0002-8758-3952), Technische Universiteit Delft

## License
The contents in the `docs/` directory are licensed under a **CC-BY 4.0** (see [CC-BY-4.0](LICENSES/CC-BY-4.0.txt) file). The source code, example notebooks, and any other file in this repository are licensed under an **Apache License v2.0** (see [Apache-License-v2.0](LICENSES/Apache-License-v2.0.txt) file).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program “LBoW”. LBoW is a python package to solve Linear Buoyancy Wave problems written by the Author(s).  
Henri Werij, Dean of Faculty of Aerospace Engineering, Technische Universiteit Delft.

&copy; 2022, D. Allaerts

## References
Linear theory of buoyancy waves is fairly standard text book material and can be found in, e.g.,
- Smith, R.B., 1980. Linear theory of stratified hydrostatic flow past an isolated mountain. Tellus 32, 348–364.
- Gossard, E., Hooke, W.H., 1975. Waves in the atmosphere, Developments in atmospheric science 2. Elsevier Scientific Publishing Company.
- Gill, A.E., 1982. Atmosphere-ocean dynamics, Nachdr. ed, International geophysics series. Acad. Press, San Diego.
- Nappo, C.J., 2002. An introduction to atmospheric gravity waves, International geophysics series. Academic Press, San Diego.
- Pedlosky, J., 2003. Waves in the ocean and atmosphere: Introduction to wave dynamics. Springer, Berlin ; New York.

For more detailed references, check out the [documentation](https://driesallaerts.github.io/lbow/intro.html).

## Would you like to contribute?

If you have any comments, feedback, or recommendations, feel free to **reach out** by sending an email to d.j.n.allaerts-1@tudelft.nl

If you would like to contribute directly, you are welcome to **fork** this repository.

Thank you and enjoy!
