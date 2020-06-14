# Helper functions
using DG_Playground, SparseArrays, LinearAlgebra
include("../examples/diffusion_utils.jl")

using Plots, DifferentialEquations, JLD2, Printf

# Mesh Stuff
K = 16     # Number of elements
n = 2      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
𝒢 = Mesh(K, n, xmin, xmax) # Generate Uniform Mesh
∇ = Gradient(𝒢) # Define a Gradient

# Define Initial Condition
u = @. exp(-2 * (xmax-xmin) / 3 * (𝒢.x - (xmax-xmin)/2)^2)

# Define hyperbolic flux
α = 0.0 # Rusanov prameter
flux_type = Rusanov(α)
field_bc = Dirichlet(0.0, 1.0)
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = 0.0 # Rusanov parameter
flux_type = Rusanov(α)
field_bc = FreeFlux()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Diffusion parameters
dt = cfl_diffusive(𝒢, 1.0) # CFL timestep
tspan  = (0.0, 5.0)
params = (∇, Φ, ∇Φ)
rhs! = diffusion!

###
affine_operator!(x,y) = rhs!(x, y, params, 0.0)
theme(:juno)
x = copy(u)
Ax = copy(u)
A, b = build_operator(affine_operator!, 𝒢)
sparse(A)
spy(A)

###
rb = reshape(b, size(𝒢.x))
plot(𝒢.x, rb)
