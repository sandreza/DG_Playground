# Defines DG data structures for convenience
# Define abstract types
include("mesh.jl")
import Base.+, Base.*, Base./ , Base.convert, Base.promote_rule, LinearAlgebra.⋅
using LinearAlgebra, SparseArrays
⋅
abstract type AbstractFlux end
abstract type AbstractGradient end
abstract type AbstractField end
abstract type AbstractFluxMethod end
abstract type AbstractBoundaryCondition end

# Structs with data
struct Gradient{𝒮} <: AbstractGradient
    grid::𝒮
end

struct Flux{𝒯, 𝒮, 𝒱} <:  AbstractFlux
    method::𝒯
    field::𝒮
    state::𝒱
end

struct Field{𝒯, 𝒮} <: AbstractField
    data::𝒯
    bc::𝒮
end

# Structs for dispatch
# Fluxes
struct NeglectFlux  <: AbstractFluxMethod end
struct Central <: AbstractFluxMethod end
struct Upwind  <: AbstractFluxMethod end

struct Rusanov{𝒯} <: AbstractFluxMethod
    α::𝒯
end

struct Slider{𝒯, 𝒮} <: AbstractFluxMethod
    α::𝒯
    v::𝒮
end

# Boundary Conditions
struct Dirichlet{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end

struct Inflow{𝒯} <: AbstractBoundaryCondition
    left::𝒯
end

struct Outflow{𝒯} <: AbstractBoundaryCondition
    left::𝒯
end

struct Neumann{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end

struct Periodic <: AbstractBoundaryCondition end
struct NoFlux   <: AbstractBoundaryCondition end

# Helper functions
function build(∇::AbstractGradient, bc::AbstractBoundaryCondition, Φ::AbstractFluxMethod; mass_matrix = false)
    #TODO build the operator in sparse representation
    return nothing
end

function compute_volume_terms(∇::AbstractArray, Φ::AbstractArray, volume_size::AbstractArray)
    q = ∇ * Φ
    @. q *= volume_size
    return q
end

function compute_volume_terms(∇::AbstractArray, Φ::AbstractField, volume_size::AbstractArray)
    q = ∇ * Φ.data
    @. q *= volume_size
    return q
end

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Central)
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

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Slider{𝒯, 𝒮}) where 𝒯 where 𝒮
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

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::AbstractBoundaryCondition, state::AbstractArray, method::NeglectFlux)
    return 𝒢.lift * zeros((𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
end

function compute_surface_terms(𝒢::AbstractMesh, Φ::AbstractField, a::Periodic, state::AbstractArray, method::Rusanov{𝒯}) where 𝒯
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
    # Now create jump in flux, (Strong-Weak form)
    @. diffs[:] -= Φ.data[𝒢.vmapM]
    # Compute Lift Operator
    lifted =  𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lifted
end

# Binary Operators
function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    # println("abstract")
    # compute volume terms
    q = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    return q
end


function ⋅(∇::AbstractGradient, Φ::Flux{𝒯, 𝒮, 𝒱}) where 𝒯 where 𝒮 where 𝒱
    # println("central")
    V = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    S = compute_surface_terms(∇.grid, Φ.field, Φ.field.bc, Φ.state, Φ.method)
    return V .+ S
end
