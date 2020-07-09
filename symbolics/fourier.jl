using Plots, FFTW
include(pwd()*"/symbolics/"*"composites.jl")

using LinearAlgebra
import Base: +, *, /, -
import LinearAlgebra: ⋅

struct FourierField{𝒯, 𝒮}
    data::𝒯
    metadata::𝒮
end

struct FourierData{𝒯}
    data::𝒯
end

-(field1::FourierData{𝒯}) where {𝒯} = -field1.data
+(field1::FourierData{𝒯}, field2::FourierData{𝒯}) where {𝒯} = field1.data + field2.data
+(field1::FourierData{𝒯}, field2::𝒮) where {𝒯, 𝒮} = field1.data + field2
+(field1::𝒯, field2::FourierData{𝒮}) where {𝒯, 𝒮} = field1 + field2.data
*(field1::FourierData{𝒯}, field2::FourierData{𝒯}) where {𝒯} = field1.data .* field2.data
*(field1::FourierData{𝒯}, field2::𝒮) where {𝒯, 𝒮} = field1.data .* field2
*(field1::𝒯, field2::FourierData{𝒮}) where {𝒯, 𝒮} = field1 .* field2.data

-(field1::FourierData{𝒯}, field2::FourierData{𝒯}) where {𝒯} = field1.data - field2.data
-(field1::FourierData{𝒯}, field2::𝒮) where {𝒯, 𝒮} = field1.data - field2
-(field1::𝒯, field2::FourierData{𝒮}) where {𝒯, 𝒮} = field1 - field2.data

function eval(Φ::FourierData{𝒯}) where {𝒯}
    return Φ
end

function eval(Φ::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮}
    return Φ.data
end

+(field1::FourierField{𝒯, 𝒮}, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(field1, field2) 
*(field1::FourierField{𝒯, 𝒮}, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Multiply(field1, field2)
-(field1::FourierField{𝒯, 𝒮}, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(field1, -field2)
-(field1::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Negative(field1)

*(op::AbstractOperation, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Multiply(op, field2)
*(op::FourierField{𝒯, 𝒮}, field2::AbstractOperation) where {𝒯, 𝒮} = Multiply(op, field2)

+(op::AbstractOperation, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(op, field2)
+(op::FourierField{𝒯, 𝒮}, field2::AbstractOperation) where {𝒯, 𝒮} = Add(op, field2)

-(op::AbstractOperation, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(op, -field2)
-(op::FourierField{𝒯, 𝒮}, field2::AbstractOperation) where {𝒯, 𝒮} = Add(op, -field2)

# Calculus

function eval(e::Gradient{𝒯, 𝒰}) where {𝒯, 𝒰}
    return fourier_derivative(eval(e.operand), e.metadata.transform, e.metadata.k)
end

⋅(∇::Derivative{𝒰}, u::𝒮) where{𝒰, 𝒮} = Gradient(u, ∇.metadata)

+(field1::AbstractOperation, field2::AbstractOperation) = Add(field1,field2)
*(field1::AbstractOperation, field2::AbstractOperation) = Multiply(field1,field2)
-(field1::AbstractOperation) = Negative(field1)
-(field1::AbstractOperation, field2::AbstractOperation) = Add(field1, Negative(field2))

# Filters fields and data
abstract type AbstractFilter end
abstract type OrszagFilter <: AbstractFilter end
abstract type NoFilter <: AbstractFilter end

struct FourierMetaData{𝒮, 𝒱, ℱ, 𝒫}
    size::𝒮
    k::𝒱
    filter::ℱ
    transform::𝒫
end
