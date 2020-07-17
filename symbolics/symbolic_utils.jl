using SymbolicUtils

abstract type AbstractAnimal end

struct GenericAnimal <: AbstractAnimal end

struct Duck{T} <: AbstractAnimal
    type::T
end

struct DuckMetaData{T} 
    sound::T
end

make_sound(a::AbstractAnimal) = println("roof")

make_sound(a::Duck) = println("quack")
make_sound(a::Duck{T}) where {T <: Float64} = println("float quack")
make_sound(a::Duck{T}) where {T <: DuckMetaData} = println("meta quack")
make_sound(a::Duck{DuckMetaData{T}}) where {T <: String} = println(a.type.sound)

dog = GenericAnimal()
int_mallard = Duck(1)
float_mallard = Duck(1.0)
meta_mallard = Duck(DuckMetaData(1))
super_duck = Duck(DuckMetaData("Super Quack!!!!!"))

make_sound(dog)
make_sound(int_mallard)
make_sound(float_mallard)
make_sound(meta_mallard)
make_sound(super_duck)

##
import Base: +, *, -, show

abstract type AbstractContext end
abstract type AbstractExpression end

struct Standard <: AbstractContext end
struct NonStandard <: AbstractContext end

struct Wrapper{T, C} <: AbstractExpression
    arg1::T
    context::C
end

Base.show(io::IO, w::Wrapper{𝒮, T}) where {𝒮, T} = println(w.arg1)

struct Add{T,S,C} <: AbstractExpression
    arg1::T
    arg2::S
    context::C
end

struct Multiply{T,S,C} <: AbstractExpression
    arg1::T
    arg2::S
    context::C
end

Wrapper(a) = Wrapper(a, Standard())
Add(a, b) = Add(a, b, Standard())
Multiply(a, b) = Multiply(a, b, Standard())

+(a::AbstractExpression, b::AbstractExpression) = Add(a, b)
*(a::AbstractExpression, b::AbstractExpression) = Multiply(a, b)
+(a::AbstractExpression, b::AbstractExpression, c::NonStandard) = Multiply(a, b)
*(a::AbstractExpression, b::AbstractExpression, c::NonStandard) = Add(a, b)

eval(e::Wrapper) = eval(e.arg1)
eval(e::Add) = eval(e.arg1) + eval(e.arg2)
eval(e::Multiply) = eval(e.arg1) * eval(e.arg2)

eval(e::Add{Wrapper{S, T}, U}) where {S, T <: NonStandard, U} = eval(e.arg1) * eval(e.arg2)

eval(Add(1,2))
eval(Add(Wrapper(1), 2))
eval(Add(Wrapper(1, NonStandard()), 2))
eval(Add(Wrapper(1), Wrapper(1, NonStandard())))

expression = Wrapper(1) + Wrapper(2)

eval(+(Wrapper(1), Wrapper(1)))
eval(+(Wrapper(1), Wrapper(1), NonStandard()))
