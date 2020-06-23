# [Boundary Conditions](@id sec:se_bc)

For the study of boundary conditions it suffices to consider a single element.


# [Boundary Conditions for the Advection Equation](@id sec:advection_bc)

First consider the case where we only have advection with wavespeed ``u`` so that the continuous equations are
```math
\begin{aligned}
\partial_t \rho + \partial_x (u\rho) &= 0
\end{aligned}
```
and the discrete equations are
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j -  \mathcal{S}_{ji} (u \rho_j) &=
    - \delta_{iN} (u\rho)^*(1)
    +
     \delta_{i0}(u\rho)^*(-1)
\end{aligned}
```
How do we enforce boundary conditions in this setting? We have only one spatial derivative, thus we must make a choice on which boundary to place boundary conditions. Since the wavespeed is positive, the flow of information is from left to right, hence it well posed to impose boundary conditions on the left end of the domain ``\rho(x=-1, t)``, but not on the right endpoint ``\rho(x = 1, t)``. (Information can't flow upstream.)

Philosophically we can take a few viewpoints on enforcing a boundary conditions such as
```math
\begin{aligned}
\rho(x=-1, t) = \rho_L
\end{aligned}
```
for some ``\rho_L \in \mathbb{R}`` (here ``L`` stands for left). Mathematically we would like to determine the flux on the boundary
```math
\begin{aligned}
    (u\rho)^*(-1,t) = F( \rho_0, \rho_L)
\end{aligned}
```
in a manner that is consistent with the boundary conditions. One option is to the think of the domain as extending beyond the left endpoint ``x  = -1`` and ``\rho`` having the value of ``\rho_L`` from the left onwards. This is equivalent to introducing a **ghost-point** ``\rho^+(-1,t) = \rho_L``. The superscript ``+`` notation is typical for ghost-points or points **exterior** to our element. The points in the interior of the element are typically denoted with a superscript ``-`` such as ``\rho^-``. Our estimate of the flux can then be taken as an average of the flux calculated with the ghost point and and the flux calculated at the interior of our domain. This is known as a **central boundary flux** as is given by the formula
```math
\begin{aligned}
(u\rho)^*(-1) = \frac{u\rho^-(-1) + u\rho^+(-1)}{2}.
\end{aligned}
```

!!! warning

    It is not always the case that one uses a central flux on the boundary.


One **option** is to choose the ghost point to be the value of the solution to the left of the domain, i.e., ``\rho^+(-1) = \rho_L``, to get
```math
\begin{aligned}
(u\rho)^*(-1) = \frac{u\rho^-(-1) + u\rho^+(-1)}{2} = \frac{u \rho_0 + u \rho_L}{2}
\end{aligned}
```
where we use a central flux to estimate the flux on the boundary. It is possible to use other types of fluxes, such as upwind or Rusanov.

Alternatively, we can state that we know what the flux should be on the boundary exactly, it should just be ``u \rho_L``, this leads to **option 2**
```math
\begin{aligned}
(u\rho)^*(-1) = u \rho_L
\end{aligned}
```
In terms of our ghost point this mean that we selected ``\rho^+(-1,t) = - \rho^-(-1,t) + 2 \rho_L``, known as the **reflection principle**. This choice of ghost point essentially converts a **central boundary flux** into an **inflow boundary flux**.

Regardless of the exact choice of ghost-point, all of them must satisfy a **consistency condition**. If ``\rho^+ = \rho^- = \text{boundary condition}``, then the equation must be consistent. For example, the following ghost point is **not consistent** for all possible ``\rho_L``
```math
\begin{aligned}
    \rho^+(-1,t) = - \rho^-(-1,t) +  \rho_L
\end{aligned}
```
since
```math
\begin{aligned}
    \rho_L = - \rho_L +  \rho_L = 0
\end{aligned}
```
which is true for ``\rho_L = 0``, but not true for any other ``\rho_L``.

We will address which of these ghost-point conditions to take shortly, but first let us discuss the right endpoint ``x = 1``. What do we do for the right endpoint? After all we cannot impose any more boundary conditions, and yet we see that we still have this ``(u\rho^*)(1)`` dangling around. We cannot just set it equal to zero, since this is effectively just introducing another boundary condition, as well as being unphysical. The only option that we have is to leave it as a free endpoint given that we cannot impose anything. We do this through the use of a **transmissive boundary flux**, in this case
```math
\begin{aligned}
    (u \rho)^*(1) = u \rho^-(x = 1) = u \rho_N
\end{aligned}
```
In terms of our ghost point, this means that we take `` \rho^+(1,t) = \rho^-(1,t) = \rho_N(t)``. This is how we let the boundary be "free". The endpoint "just is what it is".

After our brief detour on the right endpoint we are now ready to discuss the choice of fluxes on the left endpoint. In practice, to be perfectly frank, either are okay. The difference between the two boundary conditions is that one wiggles around the boundary condition (option 1) and the other enforces the boundary without wiggling around it (option 2), but both are stable with RK4. This freedom of choice, given to us by DG, is both dangerous and powerful. There are a few criteria that one can use to prune amongst all possible choices and they often revolve around either conservation arguments, physics, or energy arguments. Here "energy" does not mean kinetic or potential energy, but rather arguments that rely on [Lyapunov functions](https://en.wikipedia.org/wiki/Lyapunov_function) to bound solutions to within some region in phase space. Often times quadratic Lyapunov functions are called energy.

Let us first look at the energy argument. For constant in time ``\rho_L``, it suffices to consider the **energy argument** with ``\rho_L = 0`` since we can just shift our solution ``\rho \rightarrow \rho - \rho_L`` without changing the underlying equation. The energy argument yields the following (just multiplying through by ``\rho_i`` and utilizing the algebraic properties of the operators)
```math
\begin{aligned}
    \frac{1}{2} \partial_t (\rho_i \mathcal{M}_{ij} \rho_j) = \rho_N \left(\frac{u}{2}\rho_N - (u \rho)^*(1) \right) + \rho_0 \left( (u \rho)^*(-1) - \frac{u}{2} \rho_0 \right)
\end{aligned}
```
Our choices of numerical fluxes yields
```math
\begin{aligned}
    \frac{1}{2} \partial_t (\rho_i \mathcal{M}_{ij} \rho_j)  &= - \frac{u}{2}\rho_N^2  + \frac{u}{2} \rho_0   \rho_L
    \\
    &=  - \frac{u}{2}\rho_N^2
\end{aligned}
```
for option 1 and
```math
\begin{aligned}
    \frac{1}{2} \partial_t (\rho_i \mathcal{M}_{ij} \rho_j)  &= - \frac{u}{2}\rho_N^2  +  u \rho_0 \left( \rho_L - \frac{1}{2} \rho_0 \right)
    \\
    &= -\frac{u}{2}\left(\rho_N^2 + \rho_0^2 \right) + u \rho_0 \rho_L
    \\
    &= -\frac{u}{2}\left(\rho_N^2 + \rho_0^2 \right)
\end{aligned}
```
for option 2. From this we can see that option 2 is the more stable one since it adds more dissipation. If instead we had made the choice ``\rho^+(-1) = - \alpha \rho_0 + (1+\alpha)\rho_L``, (where ``\alpha = 1`` corresponds to the second option and ``\alpha = 0`` corresponds to the first option) we would get
```math
\begin{aligned}
    \frac{1}{2} \partial_t (\rho_i \mathcal{M}_{ij} \rho_j) &= - \frac{u}{2}\rho_N^2  +  u \rho_0 \left( \frac{1+\alpha}{2}\rho_L - \frac{\alpha}{2} \rho_0 \right)
    \\
    &= -\frac{u}{2}\left(\rho_N^2 + \alpha \rho_0^2 \right) + u \frac{1+ \alpha}{2} \rho_0 \rho_L
    \\
    &= -\frac{u}{2}\left(\rho_N^2 + \alpha \rho_0^2 \right)
\end{aligned}
```
which is okay most choices of ``\alpha \geq 0``.

In the original PDE we can go through the energy argument as well. Multiplying through by ``\rho`` and integrating over the domain, we have (with ``\rho(-1) = 0``)
```math
\begin{aligned}
\partial_t \frac{1}{2} \int_E \rho^2  &= - \frac{u}{2} \int_E \partial_x\left( \rho^2 \right)
= - \frac{u}{2} \rho^2(1).
\end{aligned}
```
This is somewhat suggestive of choosing ``\alpha = 0``. The reason why it is not exact is due to the fact that the discrete equations are different than the continuous ones, hence one cannot necessarily correspond the continuous ``\rho(1)`` with the discrete ``\rho_N``. Furthermore, if we look at the discrete conservation property of ``\rho``, i.e. we integrate over the domain we get
```math
\begin{aligned}
\partial_t \int_E \rho &= - u \rho(1)
\end{aligned}
```
on the PDE level and
```math
\begin{aligned}
 \partial_t  \int_E \rho  &= (u \rho)^*(-1) - (u \rho)^*(1) = u (1-\alpha)\frac{\rho_0}{2} - u \rho_N
\end{aligned}
```
on the discretized level. The other equation is instead suggestive of choosing ``\alpha = 1`` to make the discrete approximation look like the continuous one.

This just shows that sometimes imposing a particular kind of flux can place different principles at odds with one another. That is to say, using ``\alpha = 1`` respects the algebraic structure our conservation budget, but makes the discrete energy dissipation rate algebraically different, whereas ``\alpha = 0`` makes the discrete energy dissipation rate algebraically similar, but does not respect the conservation budget in the same way.

# Boundary Conditions for the Diffusion Equation

Let us now consider the diffusion equation,
```math
\begin{aligned}
    \partial_t \rho &= \partial_x \sigma \\
    \sigma &= \partial_x \rho
\end{aligned}
```
and its discrete form
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j &=
    - \mathcal{S}_{ji}\sigma_{j}
    + \delta_{iN} \sigma^*(1)
    - \delta_{i0} \sigma^*(-1)
    \\
    \mathcal{M}_{ij} \sigma_j &= -  \mathcal{S}_{ji} \rho_j
    + \delta_{iN}\rho^*(1)
    -
    \delta_{i0} \rho^*(-1)
\end{aligned}
```
To see how we enforce [Dirichlet boundary conditions](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition) and [Neumann Boundary conditions](https://en.wikipedia.org/wiki/Neumann_boundary_condition).

Let us first do **Neumann boundary conditions**. Because we are using a first-order formulation Neumann boundary conditions for ``\rho`` means that we are imposing Dirichlet boundary conditions on `` \sigma ``. Similar to the last section this means that the fluxes ``\rho^*(1)`` and ``\rho^*(-1)`` are free to be whatever they are on the boundary. Thus we choose **transmissive boundary fluxes** for ``\rho^*`` , i.e.
```math
\begin{aligned}
    \rho^*(1) = \rho_N \text{ and } \rho^*(-1) = \rho_0
\end{aligned}
```
For ``\sigma`` we have to determine the numerical fluxes ``\sigma^*(-1)`` and ``\sigma^*(1)``, which will depend on ``\sigma_0``, ``\sigma_N``, as well as the boundary condition at the given endpoint. For this let us suppose that we want to enforce ``\sigma(-1,t) = \sigma_L`` and ``\sigma(1, t) = \sigma_R``. Then, just like in the previous section, we have a few choices that can be made with regards to how we determine the numerical flux. Upon introduction of a ghost-point and the use of a **central boundary flux**, we take
```math
\begin{aligned}
\sigma^* = \frac{\sigma^- + \sigma^+}{2}
\end{aligned}
```
Enforcing the boundary condition ``\sigma^*(-1) = \sigma_L`` or ``\sigma^*(1) = \sigma_R`` means that, with the reflection principle, we use ghost points
```math
\begin{aligned}
    \sigma^+(-1) = - \sigma_0 + 2 \sigma_L
\end{aligned}
```
and similarly for the right endpoint. As we have done with the advection example, we can philosophize about what it means to be on the boundary and take ``\sigma^+(-1)  = \sigma_L`` and ``\sigma^+(1) = \sigma_R`` and take
```math
\begin{aligned}
    \sigma^*(1) = (\sigma_N + \sigma_R)/2
\end{aligned}
```
(similarly for the left endpoint). Or we could be more flexible and interpolate between the two with the introduction of the parameter ``\alpha``, to get
```math
\begin{aligned}
    \sigma^+ = - \alpha \sigma^- +  (1+\alpha)\sigma^{\text{bc}}.
\end{aligned}
```
We could enforce the different endpoints differently and make one choice on one boundary and a different choice on the other.

Let us go through the budget argument and the energy argument just as before. For simplicity we will choose ``\sigma_L = \sigma_R = 0``. In the budget argument, choosing ``\alpha = 1``, aka, making the flux on the boundary the flux we know that it *should be* yields a discrete scheme with exactly the same **conservation property**. To see this, first note that the discrete budget for ``\rho`` is
```math
\begin{aligned}
    \partial_t \int_{E} \rho &= \sigma^*(1) - \sigma^*(-1) = (1-\alpha_{-1})\frac{\sigma_N}{2} - (1-\alpha_{+1})\frac{\sigma_0}{2}
\end{aligned}
```
For the PDE this is
```math
\begin{aligned}
\partial_t \int_{E} \rho &= \sigma(1) - \sigma(-1) = 0
\end{aligned}
```
We can also look towards the energy argument for some guidance. In the PDE we get
```math
\begin{aligned}
    \partial_t \int_E \rho^2 &= - \int_E \sigma^2 + \rho(1) \sigma(1) - \rho(-1) \sigma(-1) \\
    &= - \int_E \sigma^2 + 0
\end{aligned}
```
The discrete equations are
```math
\begin{aligned}
    \frac{1}{2} \partial_t \rho_i \mathcal{M}_{ij} \rho_j &= - \rho_i \mathcal{S}_{ji} \sigma_j + \rho_N \sigma^*(1) - \rho_0 \sigma^*(-1) \\
    \sigma_i \mathcal{M}_{ij} \sigma_j &= - \sigma_i \mathcal{S}_{ji} \rho_j + \sigma_N \rho^*(1) - \sigma_0 \rho^*(-1)
\end{aligned}
```
We use the discrete integration by parts formula,
```math
\begin{aligned}
\rho_i \mathcal{S}_{ji} \sigma_j = \sigma_i \mathcal{S}_{ij} \rho_j =  \sigma_i (\mathcal{S}_{ij} +\mathcal{S}_{ji} - \mathcal{S}_{ji} ) \rho_j  = \rho_N \sigma_N - \rho_0 \sigma_0 - \sigma_i \mathcal{S}_{ji} \rho_j
\end{aligned}
```
to place everything in one equation. (We used the algebraic property of ``\mathcal{S}_{ij}+\mathcal{S}_{ji}`` in the derivation.) The second equation, becomes
```math
\begin{aligned}
- \rho_i \mathcal{S}_{ji} \sigma_j &=
- \sigma_i \mathcal{M}_{ij} \sigma_j  + \sigma_N (\rho^*(1) - \rho_N) - \sigma_0 (\rho^*(-1) - \rho_N)
\end{aligned}
```
Thus in the discrete equation we have
```math
\begin{aligned}
    \frac{1}{2} \partial_t \rho_i \mathcal{M}_{ij} \rho_j &= -\sigma_i \mathcal{M}_{ij} \sigma_j + \sigma_N (\rho^*(1) - \rho_N) - \sigma_0 (\rho^*(-1) - \rho_N)
    \\
    &+ \rho_N \sigma^*(1) - \rho_0 \sigma^*(-1)
\end{aligned}
```
In this case we see that ``\rho^* = \rho^-``, i.e., ``\rho^*(1,t) = \rho_N`` and ``\rho^*(-1,t) = \rho_0 `` as well as ``\sigma^* = 0`` yields a structure that is algebraically similar to the PDE. Thus in the case of the diffusion equation we see that there seems to be a ``best" choice of boundary fluxes.

For **Dirichlet boundary conditions** much is the same. Here the roles of ``\sigma`` and ``\rho`` switch places, and we "impose" transmissive flux conditions for ``\sigma`` and the different variations of boundary conditions for ``\rho``. The balance budget and energy will mimic the behavior of the PDE if we choose ``\rho^* = \rho^\text{bc}``.

# Boundary Conditions for the Advection-Diffusion Equation

The PDE is
```math
\begin{aligned}
    \partial_t \rho + \partial_x (u \rho) &= \partial_x \sigma \\
    \sigma &= \partial_x \rho
\end{aligned}
```
and the discrete form is
```math
\begin{aligned}
    \partial_t \mathcal{M}_{ij} \rho_j - \mathcal{S}_{ji} (u \rho_j) &=
    - \mathcal{S}_{ji} \sigma_j
    + \delta_{iN} \sigma^*(1)
    - \delta_{i0} \sigma^*(-1)
    \\
    &\text{ } - \delta_{iN} (u \rho)^*(1)
    +
    \delta_{i0} (u \rho)^*(-1)
    \\
    \mathcal{M}_{ij} \sigma_j &= -  \mathcal{S}_{ji}\rho_j
    + \delta_{iN} \rho^*(1)
    -
    \delta_{i0} \rho^*(-1)
\end{aligned}
```

Before discussing boundary conditions let us look at the budgets and energy argument. Integrating the PDE over the domain yields
```math
\begin{aligned}
    \partial_t \int_E \rho
    &= -u \rho(1) + u \rho(-1) +  \sigma(1) - \sigma(-1)
    \\
    \int_E \sigma &= \rho(1) - \rho(-1)
\end{aligned}
```
The discrete budget is
```math
\begin{aligned}
    \partial_t \int_E \rho
    &= -(u \rho)^*(1) + (u \rho)^*(-1) +  \sigma^*(1) - \sigma^*(-1)
    \\
    \int_E \sigma &= \rho^*(1) - \rho^*(-1)
\end{aligned}
```

The energy argument (multiplying through by ``\rho``, integrating over the domain, and performing sufficient integration by parts) for the PDE yields
```math
\begin{aligned}
    \frac{1}{2} \partial_t \rho^2 &= - \int_E \sigma^2 - \frac{u}{2}\left( \rho(1)^2 - \rho(-1)^2  \right)
    + \rho(1) \sigma(1) - \rho(-1) \sigma(-1)
    \\
    \int_E \sigma^2 &= \int_E \sigma \partial_x \rho
\end{aligned}
```
The discrete energy argument yields
```math
\begin{aligned}
    \frac{1}{2} \partial_t \rho_i \mathcal{M}_{ij} \rho_j &= - \sigma_i \mathcal{M}_{ij} \sigma_j
    \\
    &+ \rho_N \left(\frac{u}{2}\rho_N - (u \rho)^*(1) \right) - \rho_0 \left(\frac{u}{2} \rho_0   - (u \rho)^*(-1) \right)
    \\
    &+ \rho_N \sigma^*(1) - \rho_0 \sigma^*(-1)
    \\
    &+ \sigma_N ( \rho^*(1) - \rho_N) - \sigma_0 \left( \rho^*(-1) - \rho_0 \right)
    \\
    \sigma_i \mathcal{M}_{ij} \sigma_j &= - \sigma_i \mathcal{S}_{ji} \rho_j     + \sigma_N \rho^*(1)
    -
    \sigma_0 \rho^*(-1)
\end{aligned}
```

Now let us talk about Dirichlet boundary conditions on this system, and specifically ``\rho(1) = \rho(-1) = 0``. Comparing the equations, we see that to mimic the algebraic structure of the budget, we must choose ``(u\rho)^* = \rho^* = 0`` and ``\sigma^* = \sigma^-``. This is because the PDE budget is
```math
\begin{aligned}
    \partial_t \int_E \rho
    &= -0 + 0 +  \sigma(1) - \sigma(-1)
    \\
    \int_E \sigma &= 0 - 0
\end{aligned}
```
and each term in the discrete budget must be zero where the PDE budget is zero. The only terms that are non-zero are the ``\sigma`` terms, which are free to be whatever they want to be since they don't have any boundary conditions.


However, comparing the PDE energy becomes
```math
\begin{aligned}
    \frac{1}{2} \partial_t \rho^2 &= - \int_E \sigma^2
    \\
    \int_E \sigma^2 &= \int_E \sigma \partial_x \rho = - \int_E (\partial_x \sigma) \rho
\end{aligned}
```
thus we must choose,
```math
\begin{aligned}
(u\rho)^* = \frac{u \rho^- + 0}{2} \text{ , } \rho^* = 0 \text{, and } \sigma^* = \sigma^-
\end{aligned}
```
In the advective flux case where are taking a central flux between what we see in the interior to the domain and what we see exterior, whereas in the terms associated with the diffusive flux we are simply setting the value of ``\rho^*`` to what it should be on the boundary. We remind the reader that ``(u\rho)^* = \frac{u \rho^- + 0}{2}`` is an abbreviation for the simultaneous prescription ``(u\rho)^*(1) = u \rho_N/2 `` and ``(u\rho)^*(-1) = u \rho_0 /2 ``.
In order to algebraically replicate the energy budget of the PDE. One point in favor of the energy argument is that it leads to a  numerically stable algorithm whereas in the budget case one has to hope that the dissipation is enough to make things not get too out of control.

Homogeneous Neumann boundary conditions for the present problem does not make too much sense, since one is just shoving material against the right wall and not letting it escape.

# Brief Discussion on Robin Boundary Conditions / Relaxation Boundary Conditions

[Robin boundary conditions](https://en.wikipedia.org/wiki/Robin_boundary_condition) are boundary conditions of the form
```math
\begin{aligned}
    \alpha \rho(1) + \beta \sigma(1) = \gamma
\end{aligned}
```
or similarly for the left endpoint. This kind of boundary condition comes up in the CLIMA Ocean code, when one has a relaxation to a surface temperature. In that context the boundary condition is
```math
\begin{aligned}
   \left. \kappa \partial_z T \right|_\text{surface} = \lambda \left( \left. T \right|_\text{surface} - T_0 \right)
\end{aligned}
```
where ``T_0`` is some prescribed temperature field.

The difficulty in enforcing this kind of boundary condition through numerical fluxes comes from interpreting what is "transmissive" and what isn't. In the particular cases of Dirichlet or Neumann boundary conditions, it was clear that the variable without boundary conditions is the one that needs transmissive fluxes; however here the boundary condition is a linear combination of Dirichlet and Neumann.

Let us focus on potential ways to enforce this boundary condition. One way is to interpret the equation as meaning
```math
\begin{aligned}
    \alpha \rho^*(1) + \beta \sigma^*(1) = \gamma.
\end{aligned}
```
Given that robin boundary conditions involve a linear combination of ``\rho`` and ``\sigma``, this is suggestive of making the transmissive flux a linear combination of ``\rho`` and ``\sigma`` as well. Thus, we will explore options of the form
```math
\begin{aligned}
\begin{bmatrix}
\alpha & \beta \\
c & d
\end{bmatrix}
\begin{bmatrix}
\rho^* \\
\sigma^*
\end{bmatrix}
=
\begin{bmatrix}
\lambda \\
c \rho^- + d \sigma^-
\end{bmatrix}
\end{aligned}
```
where ``\alpha d - \beta c \neq 0``. Meaning that we just added an auxiliary linear equation to "enforce" transmissive boundary conditions.
The choice ``c=1``, ``d = 0``, leads to
```math
\begin{aligned}
    \rho^* &= \rho^- \\
    \sigma^* &= \frac{\gamma}{\alpha} - \frac{\beta}{\alpha} \rho^-
\end{aligned}
```
This is similar to how these boundary conditions are enforced in finite volume codes where one adds a relaxation term to the appropriate grid cell.
The choice ``c = 0``, ``d = 1``, leads to
```math
\begin{aligned}
    \rho^* &= \frac{\gamma}{\beta} - \frac{\beta}{\alpha} \sigma^- \\
    \sigma^* &= \sigma^-
\end{aligned}
```
and the choice ``c = -\beta`` and ``d = \alpha`` leads to
```math
\begin{aligned}
(\alpha^2 + \beta^2) \rho^* &= \alpha \lambda + \beta^2 \rho^- - \beta \alpha \sigma^- \\
(\alpha^2 + \beta^2) \sigma^* &= \beta \lambda - \alpha \beta \rho^- + \alpha^2 \sigma^-
\end{aligned}
```
which reduced to Dirichlet boundary conditions when ``\alpha \neq 0, \beta = 0`` and Neumann boundary conditions when ``\alpha = 0, \beta \neq 0``.

At this point it is unclear which of these choices works best, or if they all lead to the same results. Certainly one can generalize by including ghost-points. At any rate it is documented here for future reference and to open up discussion.
