# [Variational Crimes](@id sec:variational_crimes)

We will use the same notation as [Single Element](@ref sec:single_element). The other thing that makes DG different from finite volume is [aliasing](https://en.wikipedia.org/wiki/Aliasing) issues.
Consider Burgers' equation, which is the nonlinear PDE
```math
\begin{aligned}
    \partial_t \rho + \partial_x f &= 0 \\
    f &= \rho^2
\end{aligned}
```
In the Discontinuous Galerkin method all quantities are interpreted as belong to the polynomial basis set, including the flux term ``f``, i.e., we represent our functions as
```math
\begin{aligned}
    \rho &= \rho_i \ell_i \text{ and }
    f = f_i \ell_i ,
\end{aligned}
```
which are polynomials of degree ``N``, for ``i = 0, ..., N``. But yet ``f`` involves a \textit{product} of two polynomials of degree ``N``, which is a polynomial of degree ``2N``. This is not an issue in finite volume codes where ``N = 0`` so that ``2N = 0``, but can sometimes lead to issues in Discontinuous Galerkin codes. This is the classic [closure problem](https://en.wikipedia.org/wiki/Turbulence_modeling) but in the context of numerical methods.

Thus we must define in what sense we are representing a polynomial of degree ``2N`` by a polynomial of degree ``N``. We will discuss four (non-exhaustive) options
0. Exact quadrature version 1
0. Exact quadrature version 2
0. Ignoring higher order terms, (Polynomials of Degree ``N+1`` and higher)
0. Multiplying nodal points together


# Exact quadrature version 1

For exact quadrature we satisfy the nonlinear flux relation in the same sense that we are satisfying the PDE. We multiply through by test functions ``\psi`` and integrate,
```math
\begin{aligned}
   \int_\Omega \psi f = \int_\Omega \psi \rho^2.
\end{aligned}
```
We represent functions as ``f = f_j \ell_j`` and ``\rho_j \ell_j`` so that the equation becomes
```math
\begin{aligned}
    \int_\Omega \psi f_j \ell_j = \int_\Omega \psi \rho_j \ell_j \rho_k \ell_k.
\end{aligned}
```
We the take our test function space to be the same as our basis function space ``\psi = \ell_i`` to get
```math
\begin{aligned}
     \int_\Omega \ell_i  \ell_j f_j = \int_\Omega \ell_i  \ell_j  \ell_k \rho_j \rho_k.
\end{aligned}
```
The first term is just the entries of the mass matrix, i.e.,
```math
\begin{aligned}
  \mathcal{M}_{ij}  =  \int_\Omega \ell_i  \ell_j ,
\end{aligned}
```
But for quadratic nonlinearities we need to compute integrals like
```math
\begin{aligned}
   \mathcal{Q}^1_{ijk} \equiv \int_\Omega \ell_i(x) \ell_j(x) \ell_k(x)
\end{aligned}
```
which is a new object ``\mathcal{Q}^1``. Thus we represent the equation as follows
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j - \mathcal{S}_{ji} f_j &= \delta_{0i} f^*(x=-1,t) - \delta_{Ni} f^*(x=1,t)
    \\
    \mathcal{M}_{ij} f_j &= \mathcal{Q}_{ijk}^1 \rho_j \rho_k
\end{aligned}
```
# Exact quadrature version 2
This is similar to the previous one, except we directly represent the derivative term
```math
\begin{aligned}
   \mathcal{Q}^2_{ijk} \equiv \int_\Omega \ell_i'(x) \ell_j(x) \ell_k(x)
\end{aligned}
```
so that the discrete equations are
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j - \mathcal{Q}^2_{ijk} \rho_j \rho_k &= \delta_{0i} f^*(x=-1,t) - \delta_{Ni} f^*(x=1,t)
\end{aligned}
```
# Ignoring higher order terms

We carry through the multiplication of the polynomials, but just chop off all the terms corresponding to higher order terms. This can be done by ignoring higher Legendre modes (which is consistent with "Exact Quadrature version 1") or by ignoring the monomial terms (which is not consistent with "Exact Quadrature version 1").

# Multiplying nodal points together
One perspective is that we are satisfying the equation
```math
\begin{aligned}
    f_i \ell_i = (\rho_i \ell_i) (\rho_j \ell_j)
\end{aligned}
```
exactly at the nodal points. This leads to ``f_j = (\rho_j)^2``. If we do so the we are mixing discretizations.
Another perspective is that we approximate "exact quadrature version 1",  as
```math
\begin{aligned}
   [\mathcal{Q}_1]_{ijk} \equiv \int_\Omega \ell_i(x) \ell_j(x) \ell_k(x) \approx \mathcal{M}_{ij} \delta_{jk}
\end{aligned}
```
so that ``f_j = (\rho_j)^2``. This is also known as inexact quadrature or a variational crime. This approximation can lead to instability issues which are ameliorated either through dissipation, increasing resolution, or filters. When things are very resolved the higher order polynomial terms don't really play a role. This can either occur due ``h``-refinement, since things are locally linear, or ``p-``refinement for more complex reasons. The reasons are analogous to what happens with Fourier modes of analytic functions. Higher Fourier modes don't really matter since they decay exponentially and become negligible.

# Explicit Example

We now go through an explicit example of the differences for linear elements / polynomial order 1. Here the Lagrange interpolants are
```math
\begin{aligned}
    \ell_0(x) = (x-1)/ (-2) \text{ and } \ell_1(x) = (x+1) / 2
\end{aligned}
```
We can perform all integrals from the previous section
```math
\begin{aligned}
    \mathcal{Q}^1_{000} &= \mathcal{Q}^1_{111} =  1/2 \\
    \mathcal{Q}^1_{011} &= \mathcal{Q}^1_{101} = \mathcal{Q}^1_{110} = 1/6
    \\
    \mathcal{Q}^1_{100} &=
    \mathcal{Q}^1_{010} =
    \mathcal{Q}^1_{001} = 1/6
\end{aligned}
```
so that
```math
\begin{aligned}
    \begin{bmatrix}
    2/3 & 1/3 \\
    1/3 & 2/3
    \end{bmatrix}
    \begin{bmatrix}
    f_0 \\
    f_1
    \end{bmatrix}
    &=
        \begin{bmatrix}
    1/2 & 1/3 & 1/6 \\
    1/6 & 1/3 & 1/2
    \end{bmatrix}
        \begin{bmatrix}
    (\rho_0)^2 \\
    \rho_0 \rho_1 \\
    (\rho_1)^2
    \end{bmatrix}
\end{aligned}
```
which is
```math
\begin{aligned}
        \begin{bmatrix}
    f_0 \\
    f_1
    \end{bmatrix}
    &=
        \begin{bmatrix}
    5/6 & 1/3 & -1/6 \\
    -1/6 & 1/3 & 5/6
    \end{bmatrix}
        \begin{bmatrix}
    (\rho_0)^2 \\
    \rho_0 \rho_1 \\
    (\rho_1)^2
    \end{bmatrix}
\end{aligned}
```
For the other operator we have
```math
\begin{aligned}
    \mathcal{Q}^2_{000} &= [\mathcal{Q}_2]_{011} = -1/3
    \\
    \mathcal{Q}^2_{010} &= \mathcal{Q}^2_{001} = -1/6
    \\
    \mathcal{Q}^2_{1ij} &= -[\mathcal{Q}_2]_{0ij}  
\end{aligned}
```
Note that ``\mathcal{S}_{ji}\mathcal{Q}^1_{jks} = \mathcal{Q}^2_{iks}``. The two ``exact quadrature" versions coincide for linear elements.

Now let us check to see what it means to lop off higher  order terms
```math
\begin{aligned}
    4 \rho^2 &= 4 \left(\rho_0 \frac{x-1}{-2} + \rho_1 \frac{x+1}{2} \right)^2
    \\
    &= (\rho_0^2 + \rho_1^2 - 2 \rho_0 \rho_1) x^2
    + 2(\rho_0^2 - \rho_1^2)x + (\rho_0^2 + \rho_1^2 + 2 \rho_0 \rho_1)
    \\
    &= (\rho_0 - \rho_1)^2 x^2
    + (3 \rho_0^2 - \rho_1^2 + 2\rho_0 \rho_1 )\frac{x-1}{-2} + (3 \rho_1^2 - \rho_0^2 + 2 \rho_0 \rho_1)\frac{x+1}{2}
\end{aligned}
```
In this case we just neglect the term in front of ``x^2``. Thus ``f`` becomes
```math
\begin{aligned}
        \begin{bmatrix}
    f_0 \\
    f_1
    \end{bmatrix}
    &=
        \begin{bmatrix}
    3/4 & 1/2 & -1/4 \\
    -1/4 & 1/2 & 3/4
    \end{bmatrix}
        \begin{bmatrix}
    (\rho_0)^2 \\
    \rho_0 \rho_1 \\
    (\rho_1)^2
    \end{bmatrix}
\end{aligned}
```
If we had instead decided to take off the highest Legendre modes (as opposed to highest polynomial order), then we should reduce back to the exact quadrature formula (version 1). To see this first note that the first few legendre polynomials are
```math
\begin{aligned}
    P_0(x) = 1 \text{ , } P_1(x) = x  \text { , and } P_2(x) = \frac{3}{2} \left(x^2 - \frac{1}{3} \right)
\end{aligned}
```
Thus if we rewrite expressions with the highest polynomial order instead represented instead as Legendre modes via ``x^2 = x^2 - 1/3 + 1/3 = 2/3 P_2(x) + 1/3``, then we get
```math
\begin{aligned}
    4 \rho^2 &=
    (\rho_0 - \rho_1)^2 (x^2-1/3+1/3)
    \\
    &+ (3 \rho_0^2 - \rho_1^2 + 2\rho_0 \rho_1 )\frac{x-1}{-2} + (3 \rho_1^2 - \rho_0^2 + 2 \rho_0 \rho_1)\frac{x+1}{2}
    \\
    &= (\rho_0 - \rho_1)^2 \frac{2}{3} P_2(x)
    \\
    &+ \left(\frac{10}{3} \rho_0^2 - \frac{2}{3}\rho_1^2 + \frac{4}{3} \rho_0 \rho_1 \right)\frac{x-1}{-2} + \left(\frac{10}{3} \rho_1^2 - \frac{2}{3} \rho_0^2 + \frac{4}{3} \rho_0 \rho_1  \right)\frac{x+1}{2}
\end{aligned}
```
Neglecting the higher order Legendre modes (in this case just ``P_2(x)``, but in general much more). We indeed get (upon dividing by 4)
```math
\begin{aligned}
        \begin{bmatrix}
    f_0 \\
    f_1
    \end{bmatrix}
    &=
        \begin{bmatrix}
    5/6 & 1/3 & -1/6 \\
    -1/6 & 1/3 & 5/6
    \end{bmatrix}
        \begin{bmatrix}
    (\rho_0)^2 \\
    \rho_0 \rho_1 \\
    (\rho_1)^2
    \end{bmatrix}
\end{aligned}
```
as before.

And the last option is to represent ``f_0`` and ``f_1`` as
```math
\begin{aligned}
        \begin{bmatrix}
    f_0 \\
    f_1
    \end{bmatrix}
    &=
        \begin{bmatrix}
    1 & 0 & 0 \\
    0 & 0 & 1
    \end{bmatrix}
        \begin{bmatrix}
    (\rho_0)^2 \\
    \rho_0 \rho_1 \\
    (\rho_1)^2
    \end{bmatrix}
\end{aligned}
```
Again this just corresponds to multiplying things nodally, point by point. This option leads to code that is easy to write and is computationally cheap compared to the rest of the methods.
