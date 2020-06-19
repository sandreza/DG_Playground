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
\begin{align}
    \partial_t \int_E \psi \rho - \int_E  (\partial_x \psi) (u \rho) &=
    - \int_E (\partial_x \psi) \sigma
    + \int_{\partial E} \psi \sigma ^* - \int_{\partial E}  \psi (u \rho)^*
    \\
    \int_E \varphi \sigma &= - \int_{E} (\partial_x \varphi) \rho + \int_{\partial E } \varphi \rho^*
\end{align}
```
This form is taken as the definition of our partial differential equation in *weak form*.
