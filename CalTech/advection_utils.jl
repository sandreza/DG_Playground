include("../dg_utils/data_structures.jl")
# Define Flux Calculation

# Define wavespeed
const c = 2π   # speed of wave

function calculate_flux(x::AbstractArray)
    return c .* x
end

function calculate_flux(x::Number)
    return c * x
end

# Define right hand side of the differential equation
function advection!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]           # Gradient operator
    Φ = params[2]           # flux term
    Φ.state .= u            # use update state
    tmp =  ∇⋅Φ              # calculate (negative) tendency
    @. u̇ = -tmp             # correct and store it
    return nothing
end

# Determine timestep
function cfl(𝒢, c; α = c, CFL = 0.75)
    Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
    dt  = CFL * Δx / maximum([α, c])
    dt *= 0.5
    return dt
end
