"""
    Equations

Module defining critical types for formulating, manupulating,
and labeling/annotating balance laws.
"""
module Equations

include("Utilities.jl")

include("Domains.jl")

include("Core.jl")
include("Fields.jl")
include("Algebra.jl")
include("Differentiation.jl")
include("Integration.jl")

# include("PDESystems.jl")

# """
# Helper function
# """
# function ∂ₜ(Q)
#     ...
#     return Tendency(Q, args...)
# end

# """
# ∂ₜ Q
# """
# struct Tendency <: AbstractExpression
#     ...
#     ...
#     function Tendency(Q, args...)
#         ...
#         return new(Q, args...)
#     end
# end

# struct Source{ST} <: AbstractExpression
#     source_type::ST
#     ...
#     function Source(Q, args...)
#         ...
#         return new(source_type, args...)
#     end
# end

# """
# Helper functions for creating source terms
# """
# function S(q)
#     ...
#     return Source(q, ...)
# end

export BalanceLaw

"""
An abstract type describing a system of PDEs of the form:

∂ₜ Q = Σᵢ Tᵢ(Q),

where ∂ₜ Q is the `Tendency` and Σᵢ Tᵢ(Q) denotes a sum of
terms.
"""
abstract type AbstractPDESystem end

struct BalanceLaw{TT <: AbstractExpression, ET <: AbstractExpression} <: AbstractPDESystem
    lhs::TT
    rhs::ET
end
Base.:(==)(a::BalanceLaw, b::BalanceLaw) = isequal((a.lhs, a.rhs), (b.lhs, b.rhs))

"""
Allows us to write:

∂ₜ(Q) === S(q) - ∇⋅(F(q)) - ∇⋅(G(q, ∇q))

in code and immediate construct the `BalanceLaw`.

"""
Base.:(===)(lhs::AbstractExpression, rhs::AbstractExpression) = BalanceLaw(lhs, rhs)

# # Sketch of search functions for extracting specific terms
# function get_terms!(bl::BalanceLaw, terms, term_type)
#     if term_type == "Tendency"
#         return append!(terms, [bl.tendency])
#     else
#         get_terms!(bl.termsum, terms, term_type)
#     return terms
# end
# function get_terms!(expr::Operator, terms, term_type)
#     if term_type == expr.term_label
#         append!(terms, [expr])
#     end
#     for term ∈ expr.operands
#         get_terms!(term, terms, term_type)
#     end
#     return terms
# end
# # Repeat until reach Terminal nodes
# function get_terms!(expr::Terminal, terms, term_type)
#     if term_type == expr.term_label
#         append!(terms, [expr])
#     end
#     return terms
# end

# ∂ₜ q === S(q) - ∇⋅(F(q); rate=...) - ∇⋅(G(q, ∇q); rate=...)

# BoundaryCondition(q⋅n === h(x,y), "on_boundary")

"""
Sample equation:

∂ₜ q = S(q) - ∇⋅(F(q)) - ∇⋅(G(q, ∇q))                                     (eq:foo)

q - state (ρ, ρu, ρe)
F - flux of q,
G - flux of q which also depends on ∇q
S - source

When we go to DG, (eq:foo) becomes (cell-wise integral):

∫ ϕ ⋅ ∂ₜ q dx = ∫ ϕ ⋅ S(q) dx + ∫ ∇ϕ ⋅ F(q) dx - ∮ ϕ ⋅ H₁(q) ds
                + ∫ ∇ϕ ⋅ G(q) dx - ∮ ϕ ⋅ H₂(q, σ) ds,             ∀ ϕ,    (eq:DG-1)

∫ ϕ ⋅ σ dx    = -∫ ∇ϕ ⋅ g(q) dx + ∮ ϕ ⋅ H₃(g(q)) ds,              ∀ ϕ,    (eq:DG-2)

where g is some simple map (coefficient scaling) and H₃ is the numerical flux
for the auxiliary equation. (eq:DG-2) is introduced as an auxiliary variable
for approximating σ = g(∇q).
"""

# Field Signature
abstract type AbstractSignature end

struct Signature{TS, DS, RS, M} <: AbstractSignature
    time_scale::TS
    domain_space::DS
    range_space::RS
    model::M
end

# challenges
# - how to "name" subexpressions
#   - numerical fluxes
#   - boundary conditions
#   - time rates
#   - Computational performance:
#     - communication/computation (fluxes!)

end
