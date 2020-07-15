abstract type AbstractData <: AbstractExpression end
struct Data{𝒯} <: AbstractData
    data::𝒯
end

for unary_operator in unary_operators
    b_symbol = Meta.parse.(unary_operator[2]) #broadcast
    @eval -(field1::AbstractData) where {𝒯} = broadcast($b_symbol, field1.data)
end

for binary_operator in [binary_operators..., ["Negative", "-"]]
    b_symbol = Meta.parse.(binary_operator[2]) #broadcast
    @eval $b_symbol(field1::AbstractData, field2::AbstractData) = broadcast($b_symbol, field1.data, field2.data)
    @eval $b_symbol(field1::AbstractData, field2::𝒮) where {𝒮} =  broadcast($b_symbol,field1.data, field2)
    @eval $b_symbol(field1::𝒯, field2::AbstractData) where {𝒯} = broadcast($b_symbol, field1, field2.data)
end

# otherwise there is a method error
*(field1::AbstractData, field2::𝒮) where {𝒮  <: Number} = field1.data .* field2
*(field1::𝒯, field2::AbstractData) where {𝒯 <: Number} = field1 .* field2.data