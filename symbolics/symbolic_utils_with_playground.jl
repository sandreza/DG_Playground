include(pwd() * "/symbolics/abstract_core.jl")

using SymbolicUtils
import SymbolicUtils: Chain, Postwalk, Sym, Term, istree, operation, arguments, to_symbolic

SymbolicUtils.istree(a::Add) = true
SymbolicUtils.arguments(a::Add) = [a.term1, a.term2]
SymbolicUtils.operation(a::Add) = +;
SymbolicUtils.symtype(a::Add) = Number

SymbolicUtils.istree(a::Multiply) = true
SymbolicUtils.arguments(a::Multiply) = [a.term1, a.term2]
SymbolicUtils.operation(a::Multiply) = *;
SymbolicUtils.symtype(a::Multiply) = Number

a = Field(1, nothing)
c = a + a


struct Wrapper{T} <: AbstractExpression
    s::T
end

a = Wrapper(1)
c = a + a

to_expr(t::Term) = Expr(:call, operation(t), to_expr.(arguments(t))...)
to_expr(x) = x
# This is absolutely necessary
SymbolicUtils.show_simplified[] = false
symbolic_c = SymbolicUtils.to_symbolic(c);
(@rule ~b + ~b => ~b * ~b)(symbolic_c)
# to go back
to_expr(t::Term) = Expr(:call, operation(t), to_expr.(arguments(t))...)
to_expr(x) = x
eval(to_expr(symbolic_c))

# (@rule b + b => b * b)(SymbolicUtils.to_symbolic(c))
# ∂x = Sym{FnType{Tuple{Vararg{Any}}, Number}}(:∂x)
#(@rule +(~a,~b) => Add(~b, ~a))(Wrapper(1) + Wrapper(2))s