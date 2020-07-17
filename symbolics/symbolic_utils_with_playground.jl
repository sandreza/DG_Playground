include(pwd() * "/symbolics/abstract_core.jl")

using SymbolicUtils

SymbolicUtils.istree(a::Add) = true
SymbolicUtils.arguments(a::Add) = [a.term1, a.term2]
SymbolicUtils.operation(a::Add) = +;
SymbolicUtils.symtype(a) = Number

SymbolicUtils.istree(a::Multiply) = true
SymbolicUtils.arguments(a::Multiply) = [a.term1, a.term2]
SymbolicUtils.operation(a::Multiply) = *;
SymbolicUtils.symtype(a) = Number

a = Field(1, nothing)
c = a + a

struct MySymbol <: AbstractExpression
    s::Symbol
end

struct Wrapper{T} <: AbstractExpression
    s::T
end

a = Wrapper(1)
c = a + a

SymbolicUtils.istree(ex::AbstractExpression) = true
SymbolicUtils.operation(::AbstractExpression) = nothing
SymbolicUtils.arguments(::AbstractExpression) = nothing
SymbolicUtils.to_symbolic(s::AbstractExpression) = SymbolicUtils.Sym(s)
symbolic_c = SymbolicUtils.to_symbolic(c)

# (@rule b + b => b * b)(SymbolicUtils.to_symbolic(c))
# ∂x = Sym{FnType{Tuple{Vararg{Any}}, Number}}(:∂x)
