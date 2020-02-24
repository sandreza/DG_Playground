include("data_structures.jl")
include("utils.jl")
include("mesh.jl")
using Plots

# Mesh Stuff
K = 16
n = 2
xmin = 0.0
xmax = 2π
𝒢 = Mesh(K, n, xmin, xmax)

# Operator
∇ = Gradient(𝒢)

# Field
field_data = sin.(𝒢.x)
field_bc = Periodic()

# Flux
flux_type  = Central()
flux_field = Field(field_data, field_bc)
# Flux
Φ = Flux(flux_type, flux_field, field_data)


# Compute Gradient
flux_divergence = ∇⋅Φ

# Plot
theme(:juno)
plot(𝒢.x, field_data)
plot!(𝒢.x, flux_divergence, legend = false, linewidth = 3)


###
using DifferentialEquations
# Define time-stepping functions, this could be done better

# Evolve Forward in time
function solveAdvection!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]           # Gradient operator
    Φ = params[2]           # flux term
    update_flux!(Φ, u)      # use update rule
    tmp =  ∇⋅Φ              # calculate (negative) tendency
    @. u̇ = -tmp             # correct and store it
    return nothing
end

# Define time-stepping parameters

const c = 2π # speed of wave

# Determine timestep
Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
CFL = 0.75
α = 10.0 #Rusonov Flux Parameters
dt  = CFL * Δx / maximum([α, c])
dt *= 0.5 / 1

# Initial condition
u = @. exp(-4 * (𝒢.x - (xmax-xmin)/2)^2)

# Redefine Data and Flux
# Field
field_data = sin.(𝒢.x)
field_bc = Periodic()
flux_type = Central()
flux_type  = Slider(0.0, [c]) #1.0 is central, 0.0 is upwind
flux_type = Rusonov(0.0)
flux_field = Field(field_data, field_bc)
# Flux and state
v = copy(u)
Φ = Flux(flux_type, flux_field, v)
# Define rule for updating flux
function update_flux!(flux::AbstractFlux, u::AbstractArray)
    @. flux.field.data = c * u # flux calculation
    @. flux.state = u          # update state
    return nothing
end

update_flux!(Φ, u)


tspan  = (0.0, 2.0)
params = (∇, Φ)
rhs! = solveAdvection!

prob = ODEProblem(rhs!, u, tspan, params);
sol  = solve(prob, Tsit5(), dt=dt, adaptive = false);

# plotting
theme(:juno)
nt = length(sol.t)
num = 20
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
#anim = @animate
for i in indices
    plt = plot(𝒢.x, sol.u[i], xlims=(xmin, xmax), ylims = (-0.1,1.1), marker = 3,    leg = false)
    plot!(𝒢.x, sol.u[1], xlims = (xmin, xmax), ylims = (-0.1,1.1), color = "red", leg = false)
    display(plt)
    #sleep(0.25)
end

# display(anim)
