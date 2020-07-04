"""
∇⋅(F_1(q))

When we go into DG, we will need to deal with
face AND volume integrals for the DifferentialOperator:

ϕ ∇⋅(F_1(q)) * dx = -∇ϕ F_1 * dx + ϕ H_1(q) * ds
"""

export Divergence, Curl, Gradient

abstract type DifferentialOperator <: Operator end

struct Divergence{T <: AbstractExpression} <: Operator
    operand::T
end

struct Curl{T <: AbstractExpression} <: Operator
    operand::T
end

struct Gradient{T <: AbstractExpression} <: Operator
    operand::T
end

# Define operators
struct Grad end
const ∇ = Grad()
(::Grad)(operand) = Gradient(operand)
(⋅)(::Grad, operand) = Divergence(operand)
(×)(::Grad, operand) = Curl(operand)
