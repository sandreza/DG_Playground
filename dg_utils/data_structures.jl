# Defines DG data structures for convenience
# Define abstract types
import Base.+, Base.*, Base./ , Base.convert, Base.promote_rule, LinearAlgebra.⋅
using LinearAlgebra, SparseArrays
⋅
abstract type AbstractFlux end
abstract type AbstractGradient end
abstract type AbstractField end
abstract type AbstractFluxMethod end

# Structs
struct Gradient{𝒮} <: AbstractGradient
    grid::𝒮
end

struct Central <: AbstractFlux end

struct Flux{𝒯, 𝒮} <:  AbstractFlux
    method::𝒯
    field::𝒮
end

struct Field{𝒯} <: AbstractField
    values::𝒯
end

# Helper functions
function build(∇::AbstractGradient, Φ::AbstractFlux)

# Binary Operators
function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    # println("abstract")
    q = ∇.grid.D * Φ.field
    @. q *= ∇.grid.rx
    return q
end


function ⋅(∇::AbstractGradient, Φ::Flux{Central, 𝒮}) where 𝒮
    # println("central")
    q = ∇.grid.D * Φ.field
    @. q *= ∇.grid.rx
    return q
end
