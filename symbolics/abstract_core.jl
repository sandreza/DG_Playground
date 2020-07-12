include(pwd() * "/symbolics/abstract_operations.jl")
include(pwd() * "/symbolics/abstract_fields.jl")
include(pwd() * "/symbolics/abstract_data.jl")


# Include Generic Evaluation Rules
eval(a::Add{𝒯, 𝒮}) where {𝒯, 𝒮} = eval(a.term1) + eval(a.term2)
eval(a::Multiply{𝒯, 𝒮}) where {𝒯, 𝒮} = eval(a.term1) * eval(a.term2)
eval(a::Negative{𝒮}) where {𝒮} = -eval(a.term)

# Data Eval
eval(Φ::AbstractData) = Φ
# Field Eval
eval(Φ::AbstractField) = Φ.data