include("data_structures.jl")
K = 8
n = 3
xmin = 0.0
xmax = 2π
𝒢 = Mesh(K, n, xmin, xmax)
∇ = Gradient(𝒢)

a = Central
Φ = Flux(a, 𝒢.x)

g = ∇⋅Φ

theme(:juno)
scatter(𝒢.x[:], g[:], xlims = (xmin, xmax), ylims = (0.0, 2.0) )
