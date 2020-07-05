using BenchmarkTools, LinearAlgebra, SymbolicUtils
test_timing = false

abstract type AbstractOperation end
abstract type UnaryOperation{𝒯} <: AbstractOperation end
abstract type BinaryOperation{𝒯, 𝒮} <: AbstractOperation end

struct Add{𝒯, 𝒮} <: BinaryOperation{𝒯, 𝒮}
    term1::𝒯
    term2::𝒮
end

struct Multiply{𝒯, 𝒮} <: BinaryOperation{𝒯, 𝒮}
    term1::𝒯
    term2::𝒮
end

struct Dot{𝒯, 𝒮} <: BinaryOperation{𝒯, 𝒮}
    term1::𝒯
    term2::𝒮
end

function eval(e::Add{Int, Int})
    return e.term1 + e.term2
end

function eval(e::Multiply{Int, Int})
    return e.term1 * e.term2
end

function eval(e::Multiply{𝒮, 𝒯}) where {𝒮, 𝒯}
    return eval(e.term1) * eval(e.term2)
end

function eval(e::Add{𝒮, 𝒯}) where {𝒮, 𝒯}
    return eval(e.term1) + eval(e.term2)
end

# The next few are necessary for speed since eval(::Int) is relatively expensive
function eval(e::Multiply{Int, 𝒯}) where {𝒯}
    return e.term1 * eval(e.term2)
end

function eval(e::Add{Int, 𝒯}) where {𝒯}
    return e.term1 + eval(e.term2)
end

function eval(e::Multiply{𝒮, Int}) where {𝒮}
    return eval(e.term1) * e.term2
end

function eval(e::Add{𝒮, Int}) where {𝒮}
    return eval(e.term1) + e.term2
end

# Adding in the specialization to Int is actually important for speed

function eval(e::Add{Number, Number})
    return e.term1 + e.term2
end

function eval(e::Multiply{Number, Number})
    return e.term1 * e.term2
end

function eval(e::Multiply{𝒮, 𝒯}) where {𝒮, 𝒯}
    return eval(e.term1) * eval(e.term2)
end

function eval(e::Add{𝒮, 𝒯}) where {𝒮, 𝒯}
    return eval(e.term1) + eval(e.term2)
end

# The next few are necessary for speed since eval(::Number) is relatively expensive
function eval(e::Multiply{Number, 𝒯}) where {𝒯}
    return e.term1 * eval(e.term2)
end

function eval(e::Add{Number, 𝒯}) where {𝒯}
    return e.term1 + eval(e.term2)
end

function eval(e::Multiply{𝒮, Number}) where {𝒮}
    return eval(e.term1) * e.term2
end

function eval(e::Add{𝒮, Number}) where {𝒮}
    return eval(e.term1) + e.term2
end



# Test
# compare function call to expression evaluation
function more_complex(a , b, c, d, e)
    return a + ( (b+c) * (d+e))
end

if test_timing
    e = Add(3, Multiply(Add(3,3),Add(3,3)))
    @btime more_complex(3, 3, 3, 3, 3)
    @btime eval(e)
end

# Testing speed in a more computationally intense environment
# Not that eval of array is array and eval of number is a number
a = ones(1000,1000)
if test_timing
    @btime a*a
    @btime eval(Multiply(a,a))
end

# Testing Struct for Delayed evaluation
# 
import Base: +, *, /, -

struct TestingField{𝒯}
    data::𝒯
end

function eval(Φ::TestingField{𝒯}) where 𝒯
    return Φ.data
end

+(field1::TestingField{𝒯}, field2::TestingField{𝒯}) where {𝒯} = Add(field1, field2) 
*(field1::TestingField{𝒯}, field2::TestingField{𝒯}) where {𝒯} = Multiply(field1, field2)

*(op::AbstractOperation, field2::TestingField{𝒯}) where {𝒯} = Multiply(op, field2)
*(op::TestingField{𝒯}, field2::AbstractOperation) where {𝒯} = Multiply(op, field2)

+(op::AbstractOperation, field2::TestingField{𝒯}) where {𝒯} = Add(op, field2)
+(op::TestingField{𝒯}, field2::AbstractOperation) where {𝒯} = Add(op, field2)

function eval(e::Multiply{TestingField{𝒯}, TestingField{𝒯}}) where 𝒯
    return e.term1.data * e.term2.data
end

function eval(e::Add{TestingField{𝒯}, TestingField{𝒯}}) where 𝒯
    return e.term1.data + e.term2.data
end

a = ones(30,30)
b = ones(30,30)
Φ1 = TestingField(a)
Φ2 = TestingField(b)

c = a * (a + b) * a + b
Φ3 = Φ1 * (Φ1 + Φ2) * Φ1 + Φ2

if test_timing
    function checking(a,b,c; size = length(c))
        BLAS.blascopy!(size, b, 1, c, 1)
        BLAS.axpy!(1.0, a, c)
        mul!(a,b,c)
        BLAS.blascopy!(size, b, 1, c, 1)
        mul!(a,b,c)
        BLAS.blascopy!(size, b, 1, c, 1)
        BLAS.axpy!(1.0, a, b)
        return nothing 
    end
    @btime a * (a + b) * a + b
    @btime eval(Φ1 * (Φ1 + Φ2) * Φ1 + Φ2)
    @btime checking(a,b,c)
end
Φ3 = TestingField(eval(Φ3))
# for efficiency might want something like a scratch space that flips types after an operation

expr = :(Φ1 * (Φ1 + Φ2) * Φ1 + Φ2)
eval(eval(expr))
expr.args[1] = :*
eval(eval(expr))

Φ3 = (Φ1 + Φ2)

function expand(e::Add{TestingField{𝒯}, TestingField{𝒯}}) where 𝒯
    Φ1 = e.term1
    Φ2 = e.term2
    return :($Φ1 + $Φ2)
end

function expand(e::Add{𝒯, 𝒮}) where {𝒯, 𝒮}
    Φ1 = e.term1
    Φ2 = e.term2
    return :($Φ1 + $Φ2)
end

function expand(e::Multiply{𝒯, 𝒮}) where {𝒯, 𝒮}
    Φ1 = e.term1
    Φ2 = e.term2
    return :($Φ1 * $Φ2)
end


expr = expand(Φ3)
eval(expr.args[2] + expr.args[3])

Φ3 = Φ1 * (Φ1 + Φ2) * Φ1 + Φ2



three = TestingField(1) + TestingField(2)

function derivative(e::Expr)
    if e.args[1] == :sin
        argument = e.args[2]
        return :(cos($argument))
    elseif e.args[1] == :cos
        argument = e.args[2]
        return :(-sin($argument))
    elseif e.args[1] == :^
        power = e.args[3]
        return e = :($power * x^($power-1))
    end
    return nothing
end

x = 3
e = :(x^4)
derivative(e)

