##############################
# Boundary Conditions        #
##############################
function compute_boundary!(diffs, data, 𝒢, a::Periodic, calculate::Function)
    # periodic functions have no boundary
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::Inflow{𝒮}, calculate::Function) where 𝒮
    uin  = -data[𝒢.vmapI] + 2 .* calculate(bc.in)
    uout =  data[𝒢.vmapO]
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::Outflow{𝒮}, calculate::Function) where 𝒮
    uin  =  data[𝒢.vmapI]
    uout =  data[𝒢.vmapO] - 2.0 .* calculate(bc.out)
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::Dirichlet{𝒮}, calculate::Function) where 𝒮
    uin  = -data[𝒢.vmapI] + 2 .* calculate(bc.left)
    uout = -data[𝒢.vmapO] + 2 .* calculate(bc.right)
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::Dirichlet2{𝒮}, calculate::Function) where 𝒮
    uin  = calculate(bc.left)
    uout = calculate(bc.right)
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::FluxBC{𝒮}, calculate::Function) where 𝒮
    uin  = bc.left
    uout = bc.right
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
    return nothing
end

function compute_boundary!(diffs, data, 𝒢, bc::FreeFlux, calculate::Function)
    uin  = data[𝒢.vmapI]
    uout = data[𝒢.vmapO]
    diffs[𝒢.mapI]  =  @. (data[𝒢.vmapI] + uin)
    diffs[𝒢.mapO]  =  @. (data[𝒢.vmapO] + uout)
end

##############################
# Numerical Fluxes           #
##############################

# Generic Default Flux that works with Neglect flux
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::AbstractBoundaryCondition, state::AbstractArray, method::NeglectFlux, calculate::Function)
    return 𝒢.lift * zeros((𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
end

# Central
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::AbstractBoundaryCondition, state::AbstractArray, method::Central, calculate::Function)
    # compute fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] - Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Boundaries
    compute_boundary!(diffs, Φ.data, 𝒢, bc, calculate)
    # Include factor of 2 for the weak-strong form
    @. diffs *= 1.0 / 2.0
    # Compute Lift Operator
    lifted = - 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end

# Rusanov
function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, bc::AbstractBoundaryCondition, state::AbstractArray, method::Rusanov{𝒯}, calculate::Function) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] + Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    # Handle Boundaries
    compute_boundary!(diffs, Φ.data, 𝒢, bc, calculate)
    # Include factor of 2 for the weak-strong form
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusanov
    @. diffs[:] += method.α * 𝒢.normals[:] .* (state[𝒢.vmapM] - state[𝒢.vmapP]) / 2.0
    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end