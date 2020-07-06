using Plots, FFTW
include("symbolics/composites.jl")
##
function fourier_nodes(a, b, N)
    return (b-a) .* collect(0:(N-1))/N .+ a
end

function fourier_wavenumbers(a, b, N)
    up = collect(0:1:N-1)
    down = collect(-N:1:-1)
    indices = up
    indices[floor(Int, N/2):end] = down[floor(Int, N/2):end]
    wavenumbers = 2π/(b-a) .* indices
    return wavenumbers
end

function fourier_derivative(y, P, k)
    tmp = copy(y)
    dy = copy(y)
    mul!(tmp, P, y)
    @. tmp *= im * k 
    ldiv!(dy, P, tmp)
    return dy
end

N = 2^4
a,b = (0,1)
x = fourier_nodes(a, b, N)
k = fourier_wavenumbers(a, b, N)
P = plan_fft(x*(1+0im))

∂ˣ(y) = fourier_derivative(y, P, k)

y = @. sin(2π*x)*(1+0im)
z = ∂ˣ(y)
rz = real.(z)
scatter(x, rz, label = "fft derivative" )
plot!(x, 2π.*cos.(2π.*x), label = "exact" )


##
using LinearAlgebra
import Base: +, *, /, -
import LinearAlgebra: ⋅

struct FourierField{𝒯, 𝒮}
    data::𝒯
    metadata::𝒮
end


function eval(Φ::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮}
    return Φ.data
end

+(field1::FourierField{𝒯, 𝒮}, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(field1, field2) 
*(field1::FourierField{𝒯, 𝒮}, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Multiply(field1, field2)

*(op::AbstractOperation, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Multiply(op, field2)
*(op::FourierField{𝒯, 𝒮}, field2::AbstractOperation) where {𝒯, 𝒮} = Multiply(op, field2)

+(op::AbstractOperation, field2::FourierField{𝒯, 𝒮}) where {𝒯, 𝒮} = Add(op, field2)
+(op::FourierField{𝒯, 𝒮}, field2::AbstractOperation) where {𝒯, 𝒮} = Add(op, field2)

function eval(e::Multiply{FourierField{𝒯, 𝒮}, FourierField{𝒯, 𝒮}}) where {𝒯, 𝒮}
    return eval(e.term1) .* eval(e.term2)
end

function eval(e::Add{FourierField{𝒯, 𝒮}, FourierField{𝒯, 𝒮}}) where {𝒯, 𝒮}
    return eval(e.term1) .+ eval(e.term2)
end

# Calculus

struct Gradient{𝒯,𝒮} <: AbstractOperation
    operand::𝒯
    metadata::𝒮
end

function eval(e::Gradient{𝒯, 𝒰}) where {𝒯, 𝒰}
    return fourier_derivative(eval(e.operand), e.metadata.transform, e.metadata.k)
end

struct Derivative{𝒯}
    metadata::𝒯
end

⋅(∇::Derivative{𝒰}, u::𝒮) where{𝒰, 𝒮} = Gradient(u, ∇.metadata)


##
# concrete implementation
N = 2^4
a,b = (0,1)
x = fourier_nodes(a, b, N)
k = fourier_wavenumbers(a, b, N)
P = plan_fft(x*(1+0im))

∂ˣ(y) = fourier_derivative(y, P, k)

y = @. sin(2π*x)*(1+0im)
z = ∂ˣ(y)
rz = real.(z)
scatter(x, rz, label = "fft derivative" )
plot!(x, 2π.*cos.(2π.*x), label = "exact" )
metadata = (N, a, b, k, P)

abstract type AbstractFilter end
abstract type OrszagFilter <: AbstractFilter end
abstract type NoFilter <: AbstractFilter end

struct FourierMetaData{𝒮, 𝒱, ℱ, 𝒫}
    size::𝒮
    k::𝒱
    filter::ℱ
    transform::𝒫
end

fourier_meta_data = FourierMetaData(N, k, NoFilter, P)
field = FourierField(y, fourier_meta_data)
∂x = Derivative(fourier_meta_data)

eval((∂x⋅(field * field) ) + field)
eval(∂x⋅(field)) .* eval(field)

#Need to distinguish what multiply means, here we want pointwise things
function eval(e::Multiply{𝒮, 𝒯}) where {𝒮, 𝒯}
    return eval(e.term1) .* eval(e.term2)
end

# now the following will work
eval(∂x⋅(field * field) * field)
