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
