
export AbstractTensorField
export ScalarField, VectorField, MatrixField

abstract type AbstractTensorField{Rank} <: Terminal end

struct FieldBase{R}
    shape::Tuple
    fields::AbstractArray{AbstractTensorField{0}, R}
end

struct ScalarField{S, T} <: AbstractTensorField{0}
    metadata::S
    function ScalarField(; metadata = nothing)
        return new{typeof(metadata)}(metadata)
    end
end

struct VectorField{S, T} <: AbstractTensorField{1}
    metadata::S
    data::T
    function VectorField(
        fields::AbstractArray{AbstractTensorField{0}, 1};
        metadata = nothing,
        data = nothing,
    )
        data = FieldBase{1}(size(fields), fields)
        return new{typeof(metadata),
                   typeof(data)}(metadata, data)
    end
end

struct MatrixField{S, T} <: AbstractTensorField{2}
    metadata::S
    data::T
    function MatrixField(
        fields::AbstractArray{AbstractTensorField{0}, 2};
        metadata = nothing,
        data = nothing,
    )
        data = FieldBase{2}(size(fields), fields)
        return new{typeof(metadata),
                   typeof(data)}(metadata, data)
    end
end

rank(::AbstractTensorField{T}) where T = T
shape(::ScalarField) = ()
shape(o::Union{VectorField, MatrixField}) = o.data.shape
