using LinearAlgebra

import Base: +, -, *, convert, promote_rule

abstract type AbstractField <: AbstractExpression end

struct Field{𝒯, 𝒮} <: AbstractField
    data::𝒯
    metadata::𝒮
end

Field() = Field(nothing, nothing)
Field(md::AbstractMetaData) = Field(nothing, md)

# Interpret Numbers as special Fields
*(a::Number, b::AbstractExpression)  = Multiply(Field(a, nothing), b)
*(a::AbstractExpression, b::Number)  = Multiply(a, Field(b, nothing))





