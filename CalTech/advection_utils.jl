# Define Flux Calculation
function update_flux!(flux::AbstractFlux, u::AbstractArray)
    @. flux.field.data = c * u # flux calculation
    @. flux.state = u          # update state
    return nothing
end

# Define right hand side of the differential equation
function advection!(u̇, u, params, t)
    # unpack params
    ∇ = params[1]           # Gradient operator
    Φ = params[2]           # flux term
    update_flux!(Φ, u)      # use update rule
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
