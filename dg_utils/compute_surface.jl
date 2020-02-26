# Generic Default Flux that works with everything
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::AbstractBoundaryCondition, state::AbstractArray, method::NeglectFlux, calculate::Function)
    return 𝒢.lift * zeros((𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
end

################################
# Periodic Boundary Conditions #
################################
# Central
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Central, calculate::Function)
    # compute fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] - Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    @. diffs *= 1.0 / 2.0
    # Handle Periodic Boundaries
    uin  = Φ.data[𝒢.vmapO]
    uout = Φ.data[𝒢.vmapI]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] - uin) / 2
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] - uout) / 2
    # Compute Lift Operator
    lifted = - 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Slider{𝒯, 𝒮}, calculate::Function) where {𝒯, 𝒮}
    # compute fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] - Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Periodic Boundaries
    uin  = Φ.data[𝒢.vmapO]
    uout = Φ.data[𝒢.vmapI]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] - uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] - uout)
    # Adds extra part
    @. diffs = -1//2 * diffs * (𝒢.normals - (1 - method.α) * abs(method.v * 𝒢.normals)/method.v)
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* diffs)
    return lifted
end

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where 𝒯
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Periodic Boundaries
    uin  = Φ.data[𝒢.vmapO]
    uout = Φ.data[𝒢.vmapI]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0
    # Handle boundary again
    uin  = state[𝒢.vmapO]
    uout = state[𝒢.vmapI]
    diffs[𝒢.mapI]  +=  @. method.α * 𝒢.normals[𝒢.mapI] * ( state[𝒢.vmapI] - uin) / 2.0
    diffs[𝒢.mapO]  +=  @. method.α * 𝒢.normals[𝒢.mapO] * ( state[𝒢.vmapO] - uout ) / 2.0
    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end


##############################
# Inflow Boundary Conditions #
##############################

# Inflow Boundary Conditions
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::Inflow{𝒮}, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Inflow Boundary Condition
    uin  = -Φ.data[𝒢.vmapI] + 2 .* calculate(bc.in)
    uout =  Φ.data[𝒢.vmapO]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0

    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end


# Inflow Boundary Conditions
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::Inflow2{𝒮}, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Inflow Boundary Condition
    uin  =  calculate(bc.in)
    uout =  Φ.data[𝒢.vmapO]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0

    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end

##############################
# Outflow Boundary Conditions #
##############################

# Outflow Boundary Conditions
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::Outflow{𝒮}, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Outflow Boundary Condition
    uin  =  Φ.data[𝒢.vmapI]
    uout = -Φ.data[𝒢.vmapO] + 2.0 .* calculate(bc.out)
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0

    # Now create jump in flux, (Strong-Weak form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end


#################################
# Dirichlet Boundary Conditions #
#################################

# Dirichlet Boundary Conditions
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::Dirichlet{𝒮}, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Inflow Boundary Condition
    uin  = -Φ.data[𝒢.vmapI] + 2 .* calculate(bc.left)
    uout = -Φ.data[𝒢.vmapO] + 2 .* calculate(bc.right)
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0

    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end

#################################
# Free flux Boundary Conditions #
#################################

# Free flux Boundary Conditions
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::FreeFlux, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Inflow Boundary Condition
    uin  = Φ.data[𝒢.vmapI]
    uout = Φ.data[𝒢.vmapO]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] + uout)
    # Central Flux
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusonov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0

    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end
