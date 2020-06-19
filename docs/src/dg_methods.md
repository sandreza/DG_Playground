# Discontinuous Galerkin

In [CLIMA](https://github.com/CliMA/ClimateMachine.jl) the [Discontinuous Galerkin](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) method serves as our spatial discretization method. It may be thought of as a combination of spectral methods and finite volume methods. The method is a higher-order generalization of a finite volume method.


## Single Element Equation

Our goal will be to understand the Discontinuous Galerkin (DG) discretization for a single element.  We will use this to illustrate the role of boundary fluxes but also to understand differences with finite volume codes.
To better illustrate the discrete implementation of the weak formulation of the conservation equation, we consider the advection-diffusion equation
```math
\begin{aligned}
    \partial_t \rho + \partial_x \left( u \rho \right) &= \partial_x \sigma \\
    \sigma &= \partial_x \rho
\end{aligned}
```
where ``u \in \mathbb{R}``, the domain is ``x \in (-1,  1) \equiv \Omega = E``, and ``\rho(x,t)`` and/or ``\sigma(x,t)`` satisfies some prescribed boundary conditions.

### Weak Form

The weak form of the equations is obtained by multiplying through each equation by *test functions* ``\psi(x)`` , ``\varphi(x)``, integrating over the domain, and integrating by parts on the terms with derivatives, to obtain
```math
\begin{aligned}
    \partial_t \int_E \psi \rho - \int_E  (\partial_x \psi) (u \rho) &=
    - \int_E (\partial_x \psi) \sigma
    + \int_{\partial E} \psi \sigma - \int_{\partial E}  \psi (u \rho)
    \\
    \int_E \varphi \sigma &= - \int_{E} (\partial_x \varphi) \rho + \int_{\partial E } \varphi \rho
\end{aligned}
```
The terms on the boundary are interpreted as numerical fluxes, typically denoted by an asterisk as follows
```math
\begin{aligned}
    \partial_t \int_E \psi \rho - \int_E  (\partial_x \psi) (u \rho) &=
    - \int_E (\partial_x \psi) \sigma
    + \int_{\partial E} \psi \sigma^* - \int_{\partial E}  \psi (u \rho)^*
    \\
    \int_E \varphi \sigma &= - \int_{E} (\partial_x \varphi) \rho + \int_{\partial E } \varphi \rho^*
\end{aligned}
```
This form is taken as the definition of our partial differential equation in *weak form*.


### Discontinuous Galerkin Approximation
In DG we approximate the spatial structure of our functions ``\rho(x,t)`` and ``\sigma(x,t)`` by a set of linearly independent polynomials, ``\ell_i(x)`` for ``i = 0, ..., N``, within each element, so that
```math
\begin{aligned}
    \rho(x,t) = \rho_i(t) \ell_i(x) \text{ and } \sigma(x,t) = \sigma_i(t) \ell_i(x)
\end{aligned}
```
where we are using Einstein summation convention for repeated indices. To reduce notational clutter we will occasionally suppress the ``x``-dependence and ``t``-dependence in the following manner
```math
\begin{aligned}
    \rho = \rho_i \ell_i \text{ and } \sigma = \sigma_i \ell_i.
\end{aligned}
```
We have ``2(N+1)`` degrees of freedom (``N+1`` for ``\rho_i`` and ``N+1`` for ``\sigma_i``), thus we must specify ``2(N+1)`` test functions for which we are satisfying the equation. In the Galerkin method we take the test functions to be the same as the basis in which we are representing our solution, i.e., ``\psi = \ell_i(x)`` for ``i = 0,..., N`` and ``\varphi = \ell_j(x)`` for ``j = 0, ... , ``.

In index notation and with Einstein summation convection, equations become (basically just replacing ``\psi = \varphi = \ell_i`` and ``\rho = \rho_j \ell_j`` and ``\sigma = \sigma_j \ell_j``)
```math
\begin{aligned}
    \partial_t \int_E \ell_i \ell_j \rho_j - \int_E  \ell_i' \ell_j (u \rho_j) &=
    - \int_E \ell_i' \ell_j \sigma_j
    + \int_{\partial E} \ell_i \sigma^* - \int_{\partial E}  \ell_i  (u \rho)^*
    \\
    \int_E \ell_i \ell_j \sigma_j &= - \int_{E} \ell_i' \ell_j \rho_j + \int_{\partial E } \ell_i \rho^*,
\end{aligned}
```
where we have introduced the prime notation to denote a derivative, e.g., ``\ell_i'`` denotes the derivative of ``\ell_i``. We see that a few operators come up over and over again, and they are given special names. The operator whose entries are given by
```math
\begin{aligned}
    \mathcal{M}_{ij} = \int_{E} \ell_i \ell_j,
\end{aligned}
```
is known as the *mass matrix*, a name borrowed from the finite element community. The operator whose entries are given by
```math
\begin{aligned}
    \mathcal{S}_{ji} = \int_{E} \ell_i' \ell_j,
\end{aligned}
```
is known as the *stiffness matrix* a name also borrowed from the finite element community. The flipping of the indices on the entries of ``\mathcal{S}`` is purposeful and corresponds to how it is ``usually" defined.

The boundary operators become
```math
\begin{aligned}
    \int_{\partial E } \ell_i \rho^* = \ell_i(1) \rho^*(1) - \ell_i(-1) \rho^*(-1)
\end{aligned}
```
and similarly for the other terms. We are abbreviating ``\rho^*(x = 1, t)`` as ``\rho^*(1)`` and we will do so for other variables as well. The boundary terms are especially relevant since they play pivotal role in how one formulates boundary conditions.
