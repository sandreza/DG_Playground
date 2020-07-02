using GraphRecipes
using Plots
default(size=(1000, 1000))

# From http://docs.juliaplots.org/latest/graphrecipes/examples/

code2 = :( 3 + 4 * 5 )
code3 = :( sin(3) + 8+ cos(8))

navier_stokes = :(
∂ᵗ(u) = -∇⋅(u ⊗ u) + ∇⋅(ν ⊙ σ);
  ∇⋅u = 0;
    σ = ∇ ⊗ u )
plot(navier_stokes, fontsize=12, shorten=0.01, axis_buffer=0.15, nodeshape=:rect)

###
burgers = :(
∂ᵗ(u) = - ∂ˣ(u * u) + ν * ∂ˣ( σ );
    σ = ∂ˣ(u)
)
plot(burgers, fontsize=15, shorten=0.01, axis_buffer=0.15)
