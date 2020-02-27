# Why?

abstract type AbstractFruit{𝒮} end
abstract type AbstractFruitCombo end

import Base: +, *, /, -
import LinearAlgebra: ⋅

struct Smoothie{𝒯, N} <: AbstractFruit{𝒯}
    fruits::NTuple{N,AbstractFruit{𝒯}}
end

function *(🍎::AbstractFruit{𝒯}, 🍌::AbstractFruit{𝒯}) where 𝒯
    return Smoothie{𝒯, N}((🍎, 🍌))
end

# A tuple of fruits
struct Fruple{N, 𝒯} <: AbstractFruitCombo
    fruits::NTuple{N, AbstractFruit{𝒯}}
end

function +(🍎::AbstractFruit{𝒯}, 🍌::AbstractFruit{𝒯}) where 𝒯
    return Fruple{2, 𝒯}((🍎,🍌))
end

function +(🍎🍎::Fruple{N, 𝒯}, 🍌🍌::Fruple{M, 𝒯}) where {N, M, 𝒯}
    return Fruple{N + M, 𝒯}((🍎🍎.fruits..., 🍌🍌.fruits...))
end

function +(🍎::AbstractFruit{𝒯}, 🍎🍌::Fruple{M, 𝒯}) where {M, 𝒯}
    return Fruple{1 + M, 𝒯}((🍎🍌.fruits..., 🍎))
end

function +(🍎🍌::Fruple{M, 𝒯}, 🍎::AbstractFruit{𝒯}) where {M, 𝒯}
    return Fruple{1 + M , 𝒯}((🍎🍌.fruits..., 🍎))
end
###

struct Apple{𝒯, 𝒮, 𝒰, 𝒱} <: AbstractFruit{𝒯}
    tag::𝒯
    delicious::𝒮
    nutricious::𝒰
    eaten::𝒱
end

struct Banana{𝒯, 𝒮, 𝒰, 𝒱, 𝒲} <: AbstractFruit{𝒯}
    tag::𝒯
    delicious::𝒮
    nutricious::𝒰
    eaten::𝒱
    peeled::𝒲
end

# test
🍎 = Apple(true, false, true)
🍌 = Banana(true,true,true,true)
🍎🍌 = 🍎+🍌

fruit_bowl = 🍎🍌 + 🍎🍌 + 🍎🍌

fruit_bowl = fruit_bowl + 🍎
