diffusivity(q) = q < 0 ? 10.0 : 0.1

function calculate_hyperbolic_flux(x::AbstractArray)
    return x
end

function calculate_parabolic_flux(x::AbstractArray)
    return diffusivity.(x) .* x
end

function calculate_hyperbolic_flux(x::Number)
    return x
end

function calculate_parabolic_flux(x::Number)
    return diffusivity(x) * x
end



# Define right hand side of the differential equation
function convective_adjustment!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]         # Gradient operator
    Φ = params[2]         # flux term
    κ∇Φ = params[3]       # diffusive state
    Φ.state .= u          # update state
    q = ∇⊗Φ               # calculate gradient
    @. κ∇Φ.state = q      # store flux
    tmp =  ∇⋅κ∇Φ          # calculate tendency
    @. u̇ = tmp            # store it
    return nothing
end

###
using DG_Playground
using Plots, DifferentialEquations, JLD2, Printf

# Mesh Stuff
K = 5      # Number of elements
n = 4      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
𝒢 = Mesh(K, n, xmin, xmax) # Generate Mesh
∇ = Gradient(𝒢)

# Define Initial Condition
u = copy(𝒢.x)

a = 1.0
b = -1.0 / 10
bc = [a b]
neumann = true
# Define hyperbolic flux
α = 0.0 # Rusanov prameter
flux_type = Rusanov(α)
if neumann
    field_bc = FreeFlux()
else
    field_bc = Dirichlet(bc...)
end
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = 0.0 # Rusanov parameter
flux_type = Rusanov(α)
if neumann
    field_bc = Dirichlet(bc...)
else
    field_bc = FreeFlux()
end
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
κ∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Diffusion parameters
courant_number = 0.1
κ_max = 10.0
Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
dt  = courant_number * Δx^2 / κ_max
tspan  = (0.0, 10.0)
params = (∇, Φ, κ∇Φ)
rhs! = convective_adjustment!

###

# Define ODE problem
prob = ODEProblem(rhs!, u, tspan, params);
# Solve it
sol  = solve(prob, Tsit5(), dt=dt, adaptive = false);

###
# Plot it
theme(:juno)
nt = length(sol.t)
num = 30 # Number of Frames
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
anim = @animate for i in indices
    plt = plot(sol.u[i], 𝒢.x, xlims=(xmin, xmax), ylims = (-0.1,2π), marker = 3,    leg = false)
    plot!(sol.u[1], 𝒢.x, xlims = (xmin, xmax), ylims = (-0.1,2π), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    plot!(xlabel = "Temperature", ylabel = "Height")
    display(plt)
    # sleep(0.25)
end

gif(anim, pwd() * "/ca_2_scratch.gif", fps = 15)
