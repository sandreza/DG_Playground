using GraphRecipes
using Plots
theme(:default)
default(size=(1000, 1000))

# From http://docs.juliaplots.org/latest/graphrecipes/examples/

code2 = :( 3 + 4 * 5 )
code3 = :( sin(3) + 8+ cos(8))

navier_stokes = :(

∂ᵗ(u) = -( ∇⋅(u ⊗ u) )  + ∇⋅(ν ⊙ σ);
  ∇⋅u = 0;
    σ = ∇ ⊗ u

)

plot(navier_stokes, fontsize=12, shorten=0.01, axis_buffer=0.15, nodeshape=:rect)

##
burgers = :(
∂ᵗ(u) = - ∂ˣ(u * u) + ν * ∂ˣ( σ );
    σ = ∂ˣ(u)
)
p1 = plot(burgers, fontsize=15, shorten=0.01, axis_buffer=0.15, title = "primitive")

###
default(size=(1000, 1000))
# Annotated Burgers
annotated_burgers = :(
∂ᵗ(u, ::Heun) = - ∂ᵗ(*(u, u, ::Filter{S}), ::Rusanov{T}) + ν * ∂ˣ( σ , ::Central);
    σ = ∂ˣ(u, ::Central)
)
p2 = plot(annotated_burgers, fontsize=10, shorten=0.01, axis_buffer=0.15, title = "annotated")

plot(p2)

sctruct TensorField
    manifold
    tangent_space

end
###
#=
# Pseudo Code
equation =
:(
∂ᵗ(u) = - ∂ˣ(u * u) + ν * ∂ˣ( σ );
    σ = ∂ˣ(u)
)
equation = interpret(equation,
spatial = DiscontinuousGalerkin(∂ˣ),
temporal = RK4(),
state = (:u, :σ)
) # automatically annotate with central fluxes everywhere
reinterpret!(equation, 1, :(∂ˣ(u * u)), ::Rusanov(0.1)))
reinterpret!(equation, 2, :(∂ˣ(u)), ::InteriorPenalty(:u)))
Ω = Domain([-1, 1], topology = periodic)
∂Ω = ∂(Ω) # boundary (empty for peridioc, [-1, 1] for wall-bounded)
grid = Mesh(Ω, elements = h, nodes = LegendreExtrema(N))
x  = grid.x
u⁰ = exp.(-x^2)
ode_problem = InitialValueProblem(u⁰, equation, grid, Δt = 0.1, adaptive = false)
evolve(ode_problem)
=#

#=
# other thoughts
perhaps include macros like
@domain Ω = [-1, 1)
to automatically infer that it is periodic or wallbounded
For example
@domain Ω = [-1, 1] ⨂ [-1, 1) ⨂ [-1, 1)
would be wall-bounded, periodic, periodic

For boundary conditions we can start by using ghost-points and inflow boundary conditions everywhere, but we may want to have some flexibility like
reinterpret!(equation, 2, :(∂ˣ(u)), boundary_flux = ([-1], ::Inflow)))
1) If no boundary conditions are supplied we should probably interpret it as a periodic domain?
2) If wall bounded-then we should be able to set up a really small system of equations and talk about the rank of the matrix to state whether or not enough boundary conditions have been supplied.
=#


"""
Mathematical concept: tensor bundles
What is the vector space that a VectorField is a vector on at each point?
This is usually trivial (the), but in physics you come across some special cases

E.g. phase-space simulations in 6D: x, y, z, px, py, pz
    Distribution function: f(x,y,z,px,py,pz)

    grad_X(f) and grad_P(f) are both vector fields, but have different tensor bundles:

    grad_X(f) is a vector on e_x, e_y, e_z
    grad_P(f) is a vecotr on e_px, e_py, e_pz
"""
