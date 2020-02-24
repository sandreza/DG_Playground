include("../dg_utils/data_structures.jl")
include("../dg_utils/utils.jl")
include("../dg_utils/mesh.jl")
include("advection_utils.jl")
using Plots, DifferentialEquations, JLD2, Printf

# Mesh Stuff
K = 16     # Number of elements
n = 2      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
𝒢 = Mesh(K, n, xmin, xmax) # Generate Mesh
∇ = Gradient(𝒢)

# Define Initial Condition and flux type
u = @. exp(-2 * (xmax-xmin) / 3 * (𝒢.x - (xmax-xmin)/2)^2)
α = 0.0 # Rusanov prameter
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state)
# Define Advection parameters
const c = 2π   # speed of wave
dt = cfl(𝒢, c, α = α) # CFL timestep
tspan  = (0.0, 2.0)
params = (∇, Φ)
rhs! = advection!

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
    plot!(𝒢.x, sol.u[1], xlims = (xmin, xmax), ylims = (-0.1,1.1), color = "red", leg = false)
    display(plt)
    #sleep(0.25)
end

relative_error = norm(sol.u[1] .- sol.u[end]) ./ norm(sol.u[end])
relative_error_string = @sprintf("%.1e", relative_error)
println("The relative error is " * relative_error_string)
