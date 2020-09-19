using DG_Playground
include(pwd() * "/examples/diffusion_utils.jl")

using Plots, DifferentialEquations, JLD2, Printf

# Mesh Stuff
K = 6       # Number of elements
n = 4        # Polynomial Order
xmin = -3000 # left endpoint of domain
xmax = 0.0   # right endpoint of domain
h = 1000 
ΔT = 10.0
𝒢 = Mesh(K, n, xmin, xmax) # Generate Mesh
∇ = Gradient(𝒢)

# Define Initial Condition
u = @. ΔT * (exp(𝒢.x / h) - exp(xmin / h)) / (1.0 - exp(xmin / h))

# Define hyperbolic flux
α = 0.0 # Rusanov prameter
flux_type = Rusanov(α)
field_bc = FreeFlux()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = 0.0 # Rusanov parameter
flux_type = Rusanov(α)
field_bc = Dirichlet(0.0, 0.0)
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Diffusion parameters
dt = cfl_diffusive(𝒢, 1.0) # CFL timestep
dt = 3000.0
tspan  = (0.0, 5.0 * 86400)
params = (∇, Φ, ∇Φ)
rhs! = diffusion!

# Define ODE problem
prob = ODEProblem(rhs!, u, tspan, params);
# Solve it
sol  = solve(prob, Euler(), dt=dt, adaptive = false);

# Plot it
##
theme(:juno)
nt = length(sol.t)
num = 20 # Number of Frames
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
anim = @animate for i in indices
    day_label = @sprintf("%.2f ", sol.t[i] ./ 86400)
    plt = plot(sol.u[i], 𝒢.x, ylims=(xmin, xmax), xlims = (0.0,ΔT), marker = 3,    leg = false)
    plot!(sol.u[1], 𝒢.x,  color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    plot!(xlabel = "Temperature", ylabel = "Depth")
    plot!(title = "Temperature at t = " * day_label * "days")
    display(plt)
    sleep(0.1)
end
gif(anim, "ocean_simulat.gif")

##
plot(sol.u[2], 𝒢.x,  color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)