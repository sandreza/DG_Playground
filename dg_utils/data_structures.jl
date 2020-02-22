# Defines DG data structures for convenience
# Define abstract types
import Base.+, Base.*, Base./ , Base.convert, Base.promote_rule, LinearAlgebra.⋅
⋅
abstract type AbstractFlux end
abstract type AbstractGradient end
abstract type AbstractField end
abstract type AbstractFluxMethod end

abstract type Central <: AbstractFluxMethod end

struct Gradient{𝒮} <: AbstractGradient
    grid::𝒮
end

struct Flux{𝒯, 𝒮} <:  AbstractFlux
    method::𝒯
    field::𝒮
end

struct Field{𝒯} <: AbstractField
    values::𝒯
end

function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    q = ∇.grid.D * Φ.field
    @. q *= ∇.grid.rx
    return q
end

function ⋅(∇::AbstractGradient, Φ::Flux{Central, 𝒮}) where 𝒮
    q = ∇.grid.D * Φ.field
    @. q *= ∇.grid.rx
    return q
end
