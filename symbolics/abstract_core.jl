# Define Operators
import Base.show

# Unary Operators, (name, symbol)
unary_operators = []
push!(unary_operators, ["Negative", "-"])

# Binary Operators, (name, symbol)
binary_operators = []
push!(binary_operators, ["Add", "+"])
push!(binary_operators, ["Multiply", "*"])

# Define Algebraic Operators
include(pwd() * "/symbolics/abstract_operations.jl")
# Define Fields
include(pwd() * "/symbolics/abstract_fields.jl")
# Define Data
include(pwd() * "/symbolics/abstract_data.jl")


# Include Generic Evaluation Rules
for unary_operator in unary_operators
    b_name, b_symbol = Meta.parse.(unary_operator)
    @eval eval(a::$b_name{𝒮}) where {𝒮} = $b_symbol(eval(a.term))
end

for binary_operator in binary_operators
    b_name, b_symbol = Meta.parse.(binary_operator)
    @eval eval(a::$b_name{𝒮, 𝒯}) where {𝒮, 𝒯} = $b_symbol(eval(a.term1), eval(a.term2))

    @eval function Base.show(io::IO, operation::$b_name{𝒮, 𝒯}) where {𝒮, 𝒯}
        print(io, "(", operation.term1, $b_symbol , operation.term2, ")")
    end
end



# Data Eval
eval(Φ::AbstractData) = Φ
# Field Eval
eval(Φ::AbstractField) = Φ.data