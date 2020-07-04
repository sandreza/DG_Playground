"""
"""

export Sum

struct Sum <: Operator
    operands
end
Base.(:+)(t::AbstractExpression...) = Sum(t)
