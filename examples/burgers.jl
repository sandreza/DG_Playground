using DG_Playground
using Plots, DifferentialEquations, JLD2, Printf

# (1) Define the system
const κ = 0.001

function calculate_hyperbolic_flux(x)
    return κ .* x
end

function calculate_parabolic_flux(x)
    return x
end

function calculate_advective_flux(x)
    return - x .* x / 2
end

# Define right hand side of the differential equation
function burgers!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]           # Gradient operator
    Φ = params[2]           # flux term
    ∇Φ = params[3]          # diffusive state
    𝒜Φ = params[4]         # advection term
    Φ.state .= u            # update state
    𝒜Φ.state .= u           # update advective state
    q = ∇⋅Φ                 # calculate gradient
    ∇Φ.state .= q           # store gradient
    tmp =  ∇⋅∇Φ             # calculate tendency
    tmp += ∇⋅𝒜Φ             # add in advective contribution
    @. u̇ = tmp              # store it
    return nothing
end

# (2) Define Mesh
K = 20     # Number of elements
n = 2      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
𝒢 = Mesh(K, n, xmin, xmax) # Generate Uniform Mesh

# Define Gradient Object (sugar)
∇ = Gradient(𝒢)

# (3) Define Initial Condition
u = @.  exp(-2 * (xmax-xmin) / 3 * (𝒢.x - (xmax-xmin)/2)^2)

# (4) Annotate Fluxes
# Define hyperbolic flux (associated with diffusion)
α = -0.0 # Rusanov prameter
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = -0.0 # Rusanov parameter
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Advective flux
α = -0.5 # Rusanov parameter (negative)
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
𝒜Φ = Flux(flux_type, flux_field, state, calculate_advective_flux)

# (5) Define ODE Problem
Δx = 𝒢.x[2] - 𝒢.x[1]
dt = minimum([Δx^2 / κ * 0.1, abs(Δx / α)*0.3])
tspan  = (0.0, 20.0)
params = (∇, Φ, ∇Φ, 𝒜Φ)
rhs! = burgers!

# Define ODE problem
ode_problem = (rhs!, u, tspan, params)
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it
theme(:juno)
nt = length(sol.t)
num = 40 # Number of Frames
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
for i in indices
    plt = plot(𝒢.x, sol.u[i], xlims=(xmin, xmax), ylims = (-1.1,1.1), marker = 3,    leg = false)
    plot!(𝒢.x, sol.u[1], xlims = (xmin, xmax), ylims = (-1.1,1.1), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    # sleep(0.25)
end
