# Discontinuous Galerkin

In [CLIMA](https://github.com/CliMA/ClimateMachine.jl) the [Discontinuous Galerkin](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) method serves as our spatial discretization method. It may be thought of as a combination of spectral methods and finite volume methods. The method is a higher-order generalization of a finite volume method.

This section contains the following review of Discontinuos Galerkin methods
0. [Single Element](@ref sec:single_element)
0. [Boundary Conditions](@ref sec:se_bc)
0. [Variational Crimes](@ref sec:variational_crimes)
