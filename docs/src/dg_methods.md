# Discontinuous Galerkin

In [CLIMA](https://github.com/CliMA/ClimateMachine.jl) the [Discontinuous Galerkin](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) method serves as our spatial discretization method. It may be thought of as a combination of spectral methods and finite volume methods. The method is a higher-order generalization of a finite volume method.

```@contents
Pages = ["dg_methods.md"]
```


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
is known as the *stiffness matrix* a name also borrowed from the finite element community. The flipping of the indices on the entries of ``\mathcal{S}`` is purposeful and corresponds to how it is "usually" defined.

The boundary operators become
```math
\begin{aligned}
    \int_{\partial E } \ell_i \rho^* = \ell_i(1) \rho^*(1) - \ell_i(-1) \rho^*(-1)
\end{aligned}
```
and similarly for the other terms. We are abbreviating ``\rho^*(x = 1, t)`` as ``\rho^*(1)`` and we will do so for other variables as well. The boundary terms play a pivotal role in how one formulates boundary conditions as well as how one couples multiple elements together.

### Discrete Equations

Once we make choices for our functions ``\ell_i(x)`` we can write a set of discrete equations that represent the Discontinuous Galerkin scheme. We choose the ``\ell_i`` to be [Lagrange interpolants](https://en.wikipedia.org/wiki/Lagrange_polynomial) of a set of *nodal points* ``x_j`` for ``j = 0, ..., N``.  Being a Lagrange interpolant, by definition, means that ``\ell_i(x_j) = \delta_{ij}`` where ``\delta_{ij}`` is the [Kronecker delta](https://en.wikipedia.org/wiki/Kronecker_delta). The nodal points ``x_j`` are chosen as the *extrema* of [Jacobi polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials). These points are able to be efficiently calculated with either explicit formulas (as is the case with *Chebyshev polynomials* where ``x_j = \cos(\pi j / N)`` for ``j = 0, ..., N``) or by solving certain eigenvalue problems, as is the case for *Legendre polynomials*.

Regardless of the exact form, it is always the case that the endpoints are ``1`` and ``-1`` for polynomial orders bigger than one. We will use the convention that ``x_0 = -1`` and ``x_N = 1`` here, but with Chebyshev extrema one usually takes the opposite ordering. For polynomial order zero we take ``x_0 = 0`` and ``\ell_0 = 1`` in order to reduce back to a finite volume scheme. With this convention and the definition of our Lagrange interpolants, we have ``\ell_i(-1) = \delta_{iN}`` and ``\ell_i(1) = \delta_{i0}``.

With notation and conventions now established, the discrete equations are
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j - \mathcal{S}_{ji} (u \rho_j) &=
    - \mathcal{S}_{ji} \sigma_j
    + \delta_{iN} \sigma^*(1)
    - \delta_{i0} \sigma^*(-1)
    - \delta_{iN} (u \rho)^*(1)
    +  
    \delta_{i0} (u \rho)^*(-1)
    \\
    \mathcal{M}_{ij} \sigma_j &= -  \mathcal{S}_{ji}\rho_j
    + \delta_{iN} \rho^*(1)
    -
    \delta_{i0} \rho^*(-1)
\end{aligned}
```
We will explicitly calculate the mass and stiffness matrices in the following subsections.

### Explicit representations
For polynomial order zero we can work out all the Lagrange interpolants, extrema points, mass matrices, and stiffness matrices, by hand easily. Firstly, note that the extrema points for polynomial order ``N=0`` is ``x_0 = 0``. Note that here we have ``x_0 = x_N`` since ``N = 0``. The Lagrange interpolants are
```math
\begin{aligned}
    \ell_0(x) = 1
\end{aligned}
```
The mass and stiffness matrices are
```math
\begin{aligned}
    \mathcal{M}_{00} = [2] \text{ and } \mathcal{S}_{00} = [0]
\end{aligned}
```
Thus, polynomial order zero is equivalent to a finite volume scheme.

For polynomial order one we can work out all the Lagrange interpolants, extrema points, mass matrices, and stiffness matrices, by hand without too much effort. Firstly, note that the extrema points for polynomial order ``N=1`` is ``x_0 = -1`` and ``x_1 = 1``. The Lagrange interpolants are
```math
\begin{aligned}
    \ell_0(x) = \frac{x - 1}{-2} \text{ and } \ell_1(x) = \frac{x + 1}{2}
\end{aligned}
```
To calculate the entries of the mass matrix, one needs to calculate
```math
\begin{aligned}
    \mathcal{M}_{ij} = \int_{E} \ell_i \ell_j
\end{aligned}
```
i.e.
```math
\begin{aligned}
    \mathcal{M}_{00} &= \int_{-1}^1 \left(\frac{x - 1}{-2}\right)^2 dx = 2/3 \\
    \mathcal{M}_{01} &= \int_{-1}^1 \frac{x^2 - 1}{-4} dx = 1/3 .
\end{aligned}
```
By symmetry ``\mathcal{M}_{00} = \mathcal{M}_{11}`` and ``\mathcal{M}_{01} = \mathcal{M}_{10}`` so
```math
\begin{aligned}
    \mathcal{M} =
    \frac{1}{3}\begin{bmatrix}
    2 & 1 \\
    1 & 2
    \end{bmatrix}.
\end{aligned}
```
The stiffness matrix is obtained similarly, we must calculate the entries
```math
\begin{aligned}
    \mathcal{S}_{ji} = \int_{E} \ell_i' \ell_j
\end{aligned}
```
so
```math
\begin{aligned}
    \mathcal{S}^T_{00}  &= \int_{-1}^1 \frac{-1}{2}\left(\frac{x + 1}{2} \right) dx = -1/2
    \\
    \mathcal{S}^T_{01}  &= \int_{-1}^1 \frac{-1}{2}\left(\frac{x - 1}{-2} \right) dx = -1/2
    \\
    \mathcal{S}^T_{10}  &= \int_{-1}^1 \frac{1}{2}\left(\frac{x - 1}{-2} \right) dx = 1/2
    \\
    \mathcal{S}^T_{11}  &= \int_{-1}^1 \frac{1}{2}\left(\frac{x + 1}{2} \right) dx = 1/2
\end{aligned}
```
so that
```math
\begin{aligned}
    \mathcal{S}^T = \frac{1}{2} \begin{bmatrix}
        -1 & -1 \\
        1 & 1
    \end{bmatrix}.
\end{aligned}
```
Unlike polynomial order zero and finite volume schemes, all non-zero polynomial order DG discretizations have a non-zero stiffness matrix. Non-zero stiffness matrices play an important role in the determining the stability of the numerical discretization and adds extra complications that do not present themselves in the finite volume case.

Luckily this has been automated for any polynomial order, so we just display the results polynomial order two here. The extrema points are ``x_0 = -1``, ``x_1 = 0``, ``x_2 = 1``. The Lagrange interpolants are
```math
\begin{aligned}
    \ell_0(x) = x(x-1)/2
    \text{ , }
    \ell_1(x) = -(x+1)(x-1)
    \text{ , and }
    \ell_2(x) = x(x+1) / 2
\end{aligned}
```
The mass and stiffness matrices are
```math
\begin{aligned}
    \mathcal{M} = \frac{1}{15}
    \begin{bmatrix}
      4  & 2 & -1 \\
  2  & 16  &  2 \\
 -1 &  2 &  4
    \end{bmatrix}
    \text{ and }
    \mathcal{S}^T =
    \frac{1}{6}
    \begin{bmatrix}
 -3  & -4 &  1 \\
  4  & 0  & -4 \\
 -1  & 4 &  3 \\
    \end{bmatrix}
\end{aligned}
```
In DG_Playgroud we can extract these matrices as follows
```@repl
using DG_Playground
n = 3;
α = β = 0.0;
r = jacobiGL(α, β, n);
D = dmatrix(r, α, β, n);
V = vandermonde(r, α, β, n);
Mi = V * V';
M = inv(Mi)
S  = (M * D)'
```
### Algebraic Properties of Operators
We have seen some properties in the previous sections that are particular realizations of more algebraic properties of DG operators. Here we will collect \textbf{three} such algebraic properties.

The first is that
```math
\begin{aligned}
    \int_{E} \rho^2 = \rho_i \mathcal{M}_{ij} \rho_j
\end{aligned}
```
The proof is a one liner since ``\rho = \rho_i \ell_i`` and ``\int_{E}  (\ell_i \ell_j ) = \mathcal{M}_{ij}``. The calculation is as follows
```math
\begin{aligned}
  \int_{E} \rho^2  = \int_{E} (\rho_i \ell_i) (\rho_j \ell_j ) =  \rho_i \int_{E}  (\ell_i \ell_j ) \rho_j = \rho_i \mathcal{M}_{ij} \rho_j.
\end{aligned}
```
Furthermore, observe that that the mass matrix is symmetric, i.e., ``\mathcal{M}_{ij} = \mathcal{M}_{ji}``.

The next property is the \textbf{discrete integration by parts} formula
```math
\begin{aligned}
    \mathcal{S}_{ji} + \mathcal{S}_{ij} = \delta_{iN}\delta_{jN} - \delta_{i0}\delta_{j1}
\end{aligned}
```
This follows from
```math
\begin{aligned}
   \mathcal{S}_{ji} + \mathcal{S}_{ij} = \int_{E} \ell_{i}' \ell_j + \int_{E} \ell_j' \ell_i = \int_{\partial E } \ell_j \ell_i = \ell_j(1) \ell_i(1) - \ell_j(-1) \ell_i(-1)
\end{aligned}
```
via integration by parts, the definition of the [Lagrange interpolant](https://en.wikipedia.org/wiki/Lagrange_polynomial), and our convention that ``\ell_i(1) = \delta_{iN}`` and ``\ell_i(-1) = \delta_{i0}``.

This last algebraic property follows from the previous
```math
\begin{aligned}
    \rho_i \mathcal{S}_{ji} \rho_j = \frac{1}{2} \left[ (\rho_N)^2 - (\rho_0)^2 \right]
\end{aligned}
```
To see this write
```math
\begin{aligned}
    \mathcal{S}_{ji} = \frac{1}{2}\left(\mathcal{S}_{ji} + \mathcal{S}_{ij} \right) + \frac{1}{2}\left(\mathcal{S}_{ji} - \mathcal{S}_{ij} \right)
\end{aligned}
```
and use the fact that the anti-symmetric component vanishes `` \rho_i\left(\mathcal{S}_{ji} - \mathcal{S}_{ij} \right) \rho_j = 0``
