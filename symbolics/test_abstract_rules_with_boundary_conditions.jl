include(pwd()*"/symbolics" * "/dg_eval_rules.jl")
# bug eval(u * u * (u+u))
# just need to wrap things in data so that eval(field) = field.data
# and that field.data remains a closed system
# i.e. data * data should be data
##
abstract type AbstractBoundaryCondition end
struct BoundaryConditions{ℬ} <: AbstractBoundaryCondition
    bcs::ℬ
end
struct IntervalBoundaryCondition{ℬ, 𝒞} <: AbstractBoundaryCondition
    boundary::ℬ
    condition::𝒞
end
struct TransmissiveCondition <: AbstractBoundaryCondition end
struct ValueBoundaryCondition{𝒱} <: AbstractBoundaryCondition
    value::𝒱
end
## evaluate bc
for unary_operator in unary_operators
    b_name, b_symbol = Meta.parse.(unary_operator)
    @eval eval_bc(a::$b_name{𝒮}) where {𝒮} = $b_symbol(eval_bc(a.term))
end

for binary_operator in binary_operators
    b_name, b_symbol = Meta.parse.(binary_operator)
    @eval eval_bc(a::$b_name{𝒮, 𝒯}) where {𝒮, 𝒯} = $b_symbol(eval_bc(a.term1), eval_bc(a.term2))
end

eval_bc(x) = x
eval_bc(x::Field{S, T}) where {S <: Number, T} = x.data
function eval_bc(x::Field)
    tmp = [bc.condition.value for bc in u.metadata.state.bcs]
    return Data(tmp)
end

## info printing would be nice
function info(x::Field)
    println("The is a field object")
    println("The members of this struct are " , fieldnames(typeof(x)))
    println("Its name is " * x.metadata.method.name)
    println("--")
    for i in x.metadata.state.bcs
        println("It has a " , typeof(i.condition) )
        println("with value = " , i.condition.value ,  " at x = " ,  i.boundary)
        println("--")
    end
end

## Pretty printing
abstract type AbstractName end
struct Name{𝒩} <: AbstractName
    name::𝒩
end

function Base.show(io::IO, f::Field{D, T}) where {D <: Number, T}
    printstyled(io, f.data, color = 112)
end

function Base.show(io::IO, f::Field{D, DGMetaData{S, V, N}}) where {D <: Data, S, V, N <: AbstractName}
    printstyled(io, f.metadata.method.name, color = 170)
end

function Base.show(io::IO, ∇::Gradient{S, T}) where {S, T}
    printstyled(io, "∇", "(", color = 172)
    print( ∇.operand)
    printstyled( ")", color = 172)
end

##

# Domain and Boundary, fieldnames(typeof(∂Ω))
Ω  = IntervalDomain(0, 2π, periodic = false)
∂Ω = ∂(Ω)

bcL = IntervalBoundaryCondition(∂Ω.closure[1], ValueBoundaryCondition(1.0))
bcR = IntervalBoundaryCondition(∂Ω.closure[2], ValueBoundaryCondition(0.0))
bcs = BoundaryConditions((bcL, bcR))

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);

K = 8      # Number of elements
n = 1      # Polynomial Order
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
α = 0.2; # Rusanov parameter
field_md = DGMetaData(mesh, bcs, Name('u')); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
y_dg = Data(u0);
u̇ = Data(nothing);
u = Field(y_dg, field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
κ = 0.001 # Diffusivity Constant
##
# Burgers equation rhs
rhs = -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u))
pde_equation = [
    u̇ == -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u)),
]

##
# This is the resolution, one just needs to define an algebra on Data types
struct Checking{T}
    v::T
end
*(a::Checking, b::Checking) = Checking(broadcast(*, a.v , b.v))
a = randn(1000)
b = randn(1000)
@btime a .* b;
@btime Checking(a) * Checking(b);