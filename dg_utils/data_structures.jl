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

struct Flux{𝒯, 𝒮} <:  AbstractFlux
    method::𝒯
    field::𝒮
end

struct Field{𝒯, 𝒮} <: AbstractField
    data::𝒯
    bc::𝒮
end

# Structs for dispatch
# Fluxes
struct Central <: AbstractFluxMethod end
struct NoFlux  <: AbstractFluxMethod end
struct Rusonov <: AbstractFluxMethod end
struct Upwind  <: AbstractFluxMethod end


# Boundary Conditions
struct Dirichlet{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end

struct Neumann{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end

struct Periodic <: AbstractBoundaryCondition end

# Helper functions
function build(∇::AbstractGradient, bc::AbstractBoundaryCondition, Φ::AbstractFluxMethod; mass_matrix = false)
    #TODO
    return nothing
end

function compute_volume_terms(∇::AbstractArray, Φ::AbstractArray, volume_size::AbstractArray)
    q = ∇ * Φ
    @. q *= volume_size
    return q
end

function compute_volume_terms(∇::AbstractArray, Φ::Field, volume_size::AbstractArray)
    q = ∇ * Φ.data
    @. q *= volume_size
    return q
end

function compute_surface_terms(𝒢::AbstractMesh, Φ::Field, a::Periodic, method::Central)
    # compute fluxes at interface
    diffs = reshape( (Φ.data[𝒢.vmapM] - Φ.data[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    @. diffs *= 1.0 / 2.0
    # Handle Periodic Boundaries
    uin  = Φ.data[𝒢.vmapO]
    uout = Φ.data[𝒢.vmapI]
    diffs[𝒢.mapI]  =  @. (Φ.data[𝒢.vmapI] - uin) / 2
    diffs[𝒢.mapO]  =  @. (Φ.data[𝒢.vmapO] - uout) / 2
    # Compute Lift Operator
    lift = - 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* diffs)
    return lift
end

# Binary Operators
function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    # println("abstract")
    # compute volume terms
    q = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    return q
end


function ⋅(∇::AbstractGradient, Φ::Flux{𝒯, 𝒮}) where 𝒯 where 𝒮
    # println("central")
    V = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    S = compute_surface_terms(∇.grid, Φ.field, Φ.field.bc, Φ.method)
    return V .+ S
end
