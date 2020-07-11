"""
"""

export Sum

import Base: +, -

struct Sum{OT} <: Operator
    operands::OT
end

+(t::AbstractExpression...) = Sum(t)

function +(t::AbstractExpression...)
    ranks = [rank(op) for op in t]
    condition = all(ranks .== ranks[1])
    if !condition
        error("Trying to add tensors of different rank")
    end
    return Sum(t)
end

rank(a::Sum{T}) where T = rank(a.operands[1])

struct Negative{OT} <: Operator
    operands::OT
end

-(t::AbstractExpression) = Negative(t)

rank(a::Negative{T}) where T = rank(a.operand)

function -(a::AbstractExpression, b::AbstractExpression)
    return a + -b
end
