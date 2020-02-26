# Why?

abstract type AbstractFruit end
abstract type AbstractFruitCombo end

import Base: +, *, /, -
import LinearAlgebra: ⋅

expr = :(1+1)
eval(expr)
a = 3
b = 4
expr = :($a+$b)
eval(expr)


struct Apple{𝒯, 𝒮, 𝒰} <: AbstractFruit
    delicious::𝒯
    nutricious::𝒮
    eaten::𝒰
end

struct Banana{𝒯, 𝒮, 𝒰, 𝒱} <: AbstractFruit
    delicious::𝒯
    nutricious::𝒮
    eaten::𝒰
    peeled::𝒱
end

# A tuple of fruits
struct Fruple{N} <: AbstractFruitCombo where N <: Int
    fruits::NTuple{N,AbstractFruit}
end

function +(🍎::AbstractFruit, 🍌::AbstractFruit)
    return Fruple{2}((🍎,🍌))
end

function +(🍎🍎::Fruple{N}, 🍌🍌::Fruple{M}) where {N <: Int, M <: Int}
    nm = length(🍎🍎.fruits) + length(🍌🍌.fruits)
    return Fruple{nm}((🍎🍎.fruits..., 🍌🍌.fruits...))
end

function +(🍎🍎::Fruple{N}, 🍌🍌::Fruple{N}) where {N <: Number}
    nm = length(🍎🍎.fruits) + length(🍌🍌.fruits)
    return Fruple{nm}((🍎🍎.fruits..., 🍌🍌.fruits...))
end

# Needs this to find the method
function +(🍎🍎::Fruple{2}, 🍌🍌::Fruple{2})
    return Fruple{4}((🍎🍎.fruits..., 🍌🍌.fruits...))
end


###
# test
🍎 = Apple(1,2,3)
🍌 = Banana(1,2,3,4)
🍎🍌 = 🍎+🍌

🍎🍎🍌🍌 = Fruple{4}((🍎🍌.fruits..., 🍎🍌.fruits...))

nm = length(🍎🍎.fruits) + length(🍌🍌.fruits)
🍌🍌 = FruitCombo{nm}((🍎🍎.fruits..., 🍌🍌.fruits...))
🍎🍌 + 🍎🍌 + 🍎🍌
