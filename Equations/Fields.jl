
export AbstractTensorField
export ScalarField, VectorField, MatrixField

abstract type AbstractTensorField{R} <: Terminal end


struct ScalarField{S, T} <: AbstractTensorField{0}
    metadata::S
    data::T
end

struct VectorField{S, T} <: AbstractTensorField{1}
    metadata::S
    data::T
end

struct MatrixField{S, T} <: AbstractTensorField{2}
    metadata::S
    data::T
end

rank(::AbstractTensorField{T}) where T = T
