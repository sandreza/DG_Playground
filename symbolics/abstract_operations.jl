abstract type AbstractExpression end
abstract type AbstractOperation <: AbstractExpression end
abstract type UnaryOperation  <: AbstractOperation end
abstract type BinaryOperation <: AbstractOperation end

# unary operations
struct Negative{𝒯} <: UnaryOperation
    term::𝒯
end

struct PartialDerivative{𝒯, 𝒮} <: UnaryOperation
    term::𝒯
    metadata::𝒮
end

struct Add{𝒯, 𝒮} <: BinaryOperation
    term1::𝒯
    term2::𝒮
end

struct Multiply{𝒯, 𝒮} <: BinaryOperation
    term1::𝒯
    term2::𝒮
end

import Base: +, *, -

+(a::AbstractExpression, b::AbstractExpression) = Add(a, b)
*(a::AbstractExpression, b::AbstractExpression) = Multiply(a, b)
-(a::AbstractExpression) = Negative(a)
-(a::AbstractExpression, b::AbstractExpression) = Add(a, Negative(b))

