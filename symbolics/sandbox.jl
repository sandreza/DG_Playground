using SymbolicUtils, LinearAlgebra
import Base: split
import Base: +, *, /, -
import LinearAlgebra: ⋅

using SymbolicUtils: symtype, istree

# Term Labels
abstract type AbstractLabel end
abstract type Advection  <: AbstractLabel end
abstract type Diffusion  <: AbstractLabel end
abstract type AdvectiveFlux <: AbstractLabel end
abstract type DiffusiveFlux <: AbstractLabel end

# Signature
struct Signature{𝒯}
    timescale::𝒯
end
# Gradient Struct
struct Gradient{𝒟, 𝒯}
    dims::𝒟
    label::𝒯
end

# Flux Struct
struct Flux{𝒯, 𝒮, 𝒟, ℒ}
    order::𝒯
    signature::𝒮
    dims::𝒟
    label::ℒ
end

# Define some operations
function split(∇::Gradient{𝒟, 𝒯}) where {𝒟, 𝒯}
    return [Gradient(1, ∇.label[i]) for i in 1:∇.dims]
end


function ⋅(∇::Gradient, Φ::Flux)
    dims = (∇.dims, Φ.dims[1])
    return Flux(Φ.order + 1, Φ.signature,  dims, Φ.label)
end

# check some code
∇ = Gradient(3, (:x, :y, :z))
Φ = Flux(0, Signature(1), (1), Advection)
∇⋅Φ

###

struct Field{𝒮}
    space::𝒮
end

struct state{ℱ}
    fields::ℱ
end

struct Term{𝒮, 𝒪, ℒ}
    state::𝒮
    operation::𝒪
    label::ℒ
end

struct Equation{𝒯, ℒ}
    terms::𝒯
    label::ℒ
end

struct System{𝒮, ℒ}
    equations::𝒮
    label::ℒ
end
##
e = typeof(u)
@syms x::e

a = :+

3 $a 5