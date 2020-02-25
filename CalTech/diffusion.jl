include("../dg_utils/data_structures.jl")
include("../dg_utils/mesh.jl")
include("diffusion_utils.jl")

using Plots, DifferentialEquations, JLD2, Printf

# Mesh Stuff
K = 16     # Number of elements
n = 1      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
𝒢 = Mesh(K, n, xmin, xmax) # Generate Mesh
∇ = Gradient(𝒢)
∇ᴰ = Gradient(𝒢)

# Define Initial Condition
u = @. exp(-2 * (xmax-xmin) / 3 * (𝒢.x - (xmax-xmin)/2)^2)

# Define hyperbolic flux
α = 0.0 # Rusanov prameter
flux_type = Central()
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = 0.0 # Rusanov prameter
flux_type = Central()
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Diffusion parameters
dt = cfl_diffusive(𝒢, 1.0) # CFL timestep
tspan  = (0.0, 5.0)
params = (∇, Φ, ∇ᴰ, ∇Φ)
rhs! = diffusion!

# Define ODE problem
prob = ODEProblem(rhs!, u, tspan, params);
# Solve it
sol  = solve(prob, Tsit5(), dt=dt, adaptive = false);

# Plot it
theme(:juno)
nt = length(sol.t)
num = 20 # Number of Frames
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
for i in indices
    plt = plot(𝒢.x, sol.u[i], xlims=(xmin, xmax), ylims = (-0.1,1.1), marker = 3,    leg = false)
    plot!(𝒢.x, sol.u[1], xlims = (xmin, xmax), ylims = (-0.1,1.1), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    # sleep(0.25)
end
