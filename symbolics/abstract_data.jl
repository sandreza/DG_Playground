abstract type AbstractData <: AbstractExpression end

struct Data{𝒯} <: AbstractData
    data::𝒯
end

-(field1::AbstractData) where {𝒯} = -field1.data

+(field1::AbstractData, field2::AbstractData) = field1.data + field2.data
+(field1::AbstractData, field2::𝒮) where {𝒮} = field1.data + field2
+(field1::𝒯, field2::AbstractData) where {𝒯} = field1 + field2.data
*(field1::AbstractData, field2::AbstractData) = field1.data .* field2.data
*(field1::AbstractData, field2::𝒮) where {𝒮} = field1.data .* field2
*(field1::𝒯, field2::AbstractData) where {𝒯} = field1 .* field2.data

# otherwise there is a method error
*(field1::AbstractData, field2::𝒮) where {𝒮  <: Number} = field1.data .* field2
*(field1::𝒯, field2::AbstractData) where {𝒯 <: Number} = field1 .* field2.data

-(field1::AbstractData, field2::AbstractData) = field1.data - field2.data
-(field1::AbstractData, field2::𝒮) where {𝒮} = field1.data - field2
-(field1::𝒯, field2::AbstractData) where {𝒯} = field1 - field2.data