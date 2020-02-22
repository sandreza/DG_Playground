include("data_structures.jl")
include("utils.jl")
include("mesh.jl")
#using Plots
K = 8
n = 3
xmin = 0.0
xmax = 2π
𝒢 = Mesh(K, n, xmin, xmax)
∇ = Gradient(𝒢)

flux_type = Central()
# flux_type = Central
flux_field = 𝒢.x .* 𝒢.x
Φ = Flux(flux_type, flux_field)

g = ∇⋅Φ

#theme(:juno)
#scatter(𝒢.x[:], g[:], xlims = (xmin, xmax), ylims = (minimum(g), maximum(g)) )
