include("../dg_utils/data_structures.jl")

function calculate_hyperbolic_flux(x::AbstractArray)
    return x
end

function calculate_parabolic_flux(x::AbstractArray)
    return x
end

function calculate_hyperbolic_flux(x::Number)
    return x
end

function calculate_parabolic_flux(x::Number)
    return x
end

# Define right hand side of the differential equation
function diffusion!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]           # Gradient operator
    Φ = params[2]           # flux term
    ∇Φ = params[3]          # diffusive state
    Φ.state .= u            # update state
    q = ∇⊗Φ                 # calculate gradient
    ∇Φ.state .= q           # store gradient
    tmp =  ∇⋅∇Φ             # calculate tendency
    @. u̇ = tmp              # store it
    return nothing
end

# Determine timestep
function cfl_diffusive(𝒢, c; α = c, CFL = 0.1)
    Δx  = minimum(𝒢.x[2,:] -𝒢.x[1,:])
    dt  = CFL * Δx^2 / maximum([α, c])
    return dt
end
