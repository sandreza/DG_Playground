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

struct Flux{𝒯, 𝒮, 𝒱, 𝒰} <:  AbstractFlux
    method::𝒯
    field::𝒮
    state::𝒱
    calculate::𝒰
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
    in::𝒯
end

struct Outflow{𝒯} <: AbstractBoundaryCondition
    out::𝒯
end

struct Neumann{𝒯} <: AbstractBoundaryCondition
    left::𝒯
    right::𝒯
end

struct Periodic <: AbstractBoundaryCondition end
struct NoFlux   <: AbstractBoundaryCondition end
struct FreeFlux <: AbstractBoundaryCondition end

# Helper functions
function build(∇::AbstractGradient, bc::AbstractBoundaryCondition, flux_type::AbstractFluxMethod, state::AbstractArray; mass_matrix = false)
    eⁿ = zeros(length(state))
    A = eⁿ * eⁿ'
    field_data = copy(state)
    flux_field = Field(field_data, bc)
    state = copy(u)
    Φ = Flux(flux_type, flux_field, state, x->x)
    for i in eachindex(u)
        eⁿ[i] = 1.0
        state[:] .=
        Φ.field.data[:] .= eⁿ
        Φ.state[:] .= eⁿ
        tmp = ∇⋅Φ
        eⁿ[i] = 0.0
        if mass_matrix
            tmp =  𝒢.M * tmp
            @. tmp /= 𝒢.rx
            @. A[:,i] = tmp[:]
        else
            @. A[:,i] = tmp[:]
        end
    end
    return A
end

# Binary Operators
function ⋅(∇::AbstractGradient, Φ::AbstractFlux)
    q = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)
    return q
end


function ⋅(∇::AbstractGradient, Φ::Flux{𝒯, 𝒮, 𝒱, 𝒰}) where {𝒯, 𝒮, 𝒱, 𝒰}
    # calculate flux
    tmp = Φ.calculate(Φ.state)
    Φ.field.data .= tmp

    # volume terms
    V = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)

    # surface terms
    S = compute_surface_terms(∇.grid, Φ.field, Φ.field.bc, Φ.state, Φ.method, Φ.calculate)
    return V .+ S
end

# (same in 1D)
function ⊗(∇::AbstractGradient, Φ::Flux{𝒯, 𝒮, 𝒱, 𝒰}) where {𝒯, 𝒮, 𝒱, 𝒰}
    # calculate flux
    tmp = Φ.calculate(Φ.state)
    Φ.field.data .= tmp

    # volume terms
    V = compute_volume_terms(∇.grid.D, Φ.field, ∇.grid.rx)

    # surface terms
    S = compute_surface_terms(∇.grid, Φ.field, Φ.field.bc, Φ.state, Φ.method, Φ.calculate)
    return V .+ S
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

include("compute_surface.jl")
