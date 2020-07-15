
abstract type AbstractEquation end
abstract type AbstractSystem end

struct PDEEquation{TT <: AbstractExpression, ET <: AbstractExpression} <: AbstractEquation
    lhs::TT
    rhs::ET
end
Base.:(==)(a::AbstractExpression, b::AbstractExpression) = PDEEquation(a, b)

struct PDESystem <: AbstractSystem
    equations
    domain
    bcs
    initial_conditions
    metadata
end
