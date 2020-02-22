# Defines DG data structures for convenience
# Define abstract types
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

struct Field{𝒯, 𝒮, 𝒰} <: AbstractField
    data::𝒯
    bc::𝒮
    bc_type::𝒰
end

# Structs for dispatch
# Fluxes
struct Central <: AbstractFluxMethod end

# Boundary Conditions
struct Dirichlet <: AbstractBoundaryCondition end
struct Neumann <: AbstractBoundaryCondition end

# Helper functions
function build(∇::AbstractGradient, Φ::AbstractFluxMethod)
    return nothing
end

function compute_volume_terms(∇::AbstractArray, Φ::AbstractArray, volume_size::AbstractArray)
    q = ∇ * Φ
    @. q *= volume_size
    return q
end

function compute_surface_terms()
end

# Binary Operators
function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    # println("abstract")
    # compute volume terms
    q = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    compute_surface_terms()
    return q
end


function ⋅(∇::AbstractGradient, Φ::Flux{Central, 𝒮}) where 𝒮
    # println("central")
    q = ∇.grid.D * Φ.field
    @. q *= ∇.grid.rx
    return q
end
