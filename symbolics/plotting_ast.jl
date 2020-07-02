using GraphRecipes
using Plots
default(size=(1000, 1000))

# From http://docs.juliaplots.org/latest/graphrecipes/examples/

code2 = :(3 + 4 * 5)
code3 = :( sin(3) + 8+ cos(8))

expr = :( ∂ᵗ(u) = -∇⋅(u⊗u) + ∇⋅(ν ⊙ ∇u);
            ∇⋅u = 0 )
plot(expr, fontsize=12, shorten=0.01, axis_buffer=0.15, nodeshape=:rect)
