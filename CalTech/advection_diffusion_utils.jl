include("../dg_utils/data_structures.jl")
const c = 2π
const κ = 0.4

function calculate_hyperbolic_flux(x::AbstractArray)
    return κ .* x
end

function calculate_parabolic_flux(x::AbstractArray)
    return κ .* x
end

function calculate_hyperbolic_flux(x::Number)
    return x
end

function calculate_parabolic_flux(x::Number)
    return x
end

function calculate_advective_flux(x::AbstractArray)
    return -c .* x
end

function calculate_advective_flux(x::Number)
    return -c * x
end

# Define right hand side of the differential equation
function advection_diffusion!(u̇, u, params, t)
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

# Determine timestep
function cfl_diffusive(𝒢, c; α = c, CFL = 0.1)
    Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
    dt  = CFL * Δx^2 / maximum([α, c])
    return dt
end

# Determine timestep
function cfl(𝒢, c; α = c, CFL = 0.75)
    Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
    dt  = CFL * Δx / maximum([α, c])
    dt *= 0.5
    return dt
end

# Determine timestep

function cfl_advection_diffusion(𝒢, c; α = c, CFL = 0.1)
    return minimum([cfl(𝒢, c; α = α, CFL = CFL), cfl_diffusive(𝒢, c; α = α, CFL = CFL)])
end
