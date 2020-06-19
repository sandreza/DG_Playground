# Convective Adjustment

Convective adjustment is a simple parameterization that attempts the capture the effect of mixing due to convection. Physically this occurs because dense water parcels tend to sink and and light water parcels tend to rise

## Mathematical Form

Typically the effect of convective adjustment is captured via a nonlinear diffusivity such as
```math
\begin{aligned}
\kappa(\rho) = \begin{cases}
\kappa_1 & \text{ if }\text{ }  \partial_z \rho < 0
\\
\kappa_2 & \text{ otherwise}
\end{cases}
\end{aligned}
```
where `` \kappa_1 \gg \kappa_2 ``, and `` z `` is aligned with the direction of gravity. Thinking of $$ \rho $$ as density, a simple parameterization of convection is of the form
```math
\begin{aligned}
\partial_t \rho &= \nabla \cdot \left[ \kappa(\rho) \nabla \rho \right]
\end{aligned}
```

## Typical Time-Discretization

A typical time-discretization would be
```math
\begin{aligned}
\rho^{n+1} - \Delta t \partial_z \left[ \kappa(\rho^{n}) \partial_z \rho^{n+1} \right] &= \rho^{n} + \Delta t \left( f^n + \nabla^H \cdot \left[ \kappa(\rho^n) \nabla^H \rho^n \right] \right)
\end{aligned}
```
where the forcing function `` f^n `` comes from boundary condition and we have split the gradient operator into the vertically aligned component ``z`` and the other (horizontal) directions.

## Simplification

There are a wide variety of functional forms that $\kappa(\rho^n)$ can take on, but typically it is similar to
```math
\begin{aligned}
    \kappa(\rho^n) \approx
    \begin{cases}
    \kappa_1 & \text{ if } z > h \\
    \kappa_2 & \text{ if } z \leq h
    \end{cases}
\end{aligned}
```
where ``z \in [0, L]`` and ``h`` can take on all values between ``[0, L]``. The reason why there is usually just one place that ``\kappa`` changes values of diffusivity has to do typical physical scenarios that arise in the ocean / atmosphere. The ocean interior is stably stratified as there is no appreciably strong convection. In the ocean cooling comes from the surface of the ocean and leads to mixing that starts in the upper ocean and progresses towards the ocean abyss. The solution to linear systems of this form (when ``\Delta t`` becomes large), is essentially constant in the region of high diffusivity.
