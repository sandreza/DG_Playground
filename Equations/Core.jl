
export AbstractExpression
export Terminal
export Operator
export shape, operands, reconstruct
export geometric_dimension, topological_dimension
export compute_hash, hash_behavior

"""
Base type for all PDE expressions
"""
abstract type AbstractExpression end

"""
An expression that does not depend on any other expression.

Why? Expressions (PDEs) can be represented as a syntax tree
and it will be beneficial for us to explicitly define Terminal
expressions so tree visitors (functions traversing the AST)
know when they reach the end of a branch.
"""
abstract type Terminal <: AbstractExpression end

#= Do we need to distinguish between the two types of
quantities?

# Different types of `Terminal` quantities
# PrognosticQuantity like the state is a terminal quantity.
# What other things could be terminal quantities?

"""
Q = (Momentum, density, total energy, etc.)
"""
abstract type PrognosticQuantity <: Terminal end

"""
pressure/ Exner function, potential temp / temp.
vorticity, PV, etc.
"""
abstract type DiagnosticQuantity <: Terminal end
=#

"""
An expression obtained after applying an operator to
an existing expression. For example, differentiation.
"""
abstract type Operator <: AbstractExpression end

function Base.repr(o::Operator)
    str_ops = join((repr(op) for op in operands(o)), ",")
    return "$(typeof(o))$(str_ops)"
end

const VarTuple{T} = NTuple{N, T} where N
const Dimension = Int32
const DimensionTuple = VarTuple{Dimension}

function shape end 
function operands end
function reconstruct end

function reconstruct(t::Terminal, operands::VarTuple{AbstractExpression}...)
    !isempty(operands) && error("Terminal has no operands")
    return t
end 

function reconstruct(o::Operator, operands::VarTuple{AbstractExpression}...)
    return typeof(o)(operands...)
end

geometric_dimension(::AbstractExpression) = -1
topological_dimension(x::Any)::Dimension = x.topological_dimension

hash_behavior(x::Any) = x
hash_behavior(x::AbstractExpression) = x.hash_code
compute_hash(o::Operator) = hash((typecode(o), (hash(op) for op in operands(o))...))
