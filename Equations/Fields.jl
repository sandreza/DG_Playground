
export AbstractTensorField
export ScalarField #, VectorField, MatrixField
abstract type AbstractTensorField{Rank} <: Terminal end

# struct FieldBase
#     shape::Tuple
#     fields::AbstractArray{AbstractTensorField{0}}
# end

struct ScalarField{S} <: AbstractTensorField{0}
    metadata::S
end

function ScalarField(; metadata = nothing)
    return ScalarField{typeof(metadata)}(metadata)
end

# struct VectorField{S, T} <: AbstractTensorField{1}
#     metadata::S
#     data::T
# end

# function VectorField(
#     fields::T;
#     metadata = nothing,
#     data = nothing,
# ) where T
#     data = FieldBase{1}(size(fields), fields)
#     return VectorField{typeof(metadata),
#                        typeof(data)}(metadata, data)
# end

# struct MatrixField{S, T} <: AbstractTensorField{2}
#     metadata::S
#     data::T
# end

# function MatrixField(
#     fields::AbstractArray{AbstractTensorField{0}, 2};
#     metadata = nothing,
#     data = nothing,
# )
#     data = FieldBase{2}(size(fields), fields)
#     return MatrixField{typeof(metadata),
#                        typeof(data)}(metadata, data)
# end

rank(::AbstractTensorField{T}) where T = T
shape(::ScalarField) = ()
# shape(o::Union{VectorField, MatrixField}) = o.data.shape
