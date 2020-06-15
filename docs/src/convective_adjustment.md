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
where `` \kappa_2 \gg \kappa_1 ``, and `` z `` is aligned with the direciton of gravity. Thinking of $$ \rho $$ as density, a simple parameterization of convection is of the form
```math
\begin{aligned}
\partial_t \rho &= \nabla \cdot \left[ \kappa(\rho) \nabla \rho \right]
\end{aligned}
```

## Typical Time-Discretization

A typical time-discretization would be
```math
\begin{aligned}
\rho^{n+1} - \Delta t \nabla \cdot \left[ \kappa(\rho^{n}) \nabla \rho^{n+1} \right] &= 0.0
\end{aligned}
```
