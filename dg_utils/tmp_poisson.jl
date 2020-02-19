include("field.jl")

using Plots
using BenchmarkTools
using DifferentialEquations
using BandedMatrices

# choose eqn type
periodic = false #need to keep as false
timings = false   #to see how different linear solvers perform

# set number of DG elements and polynomial order
K = 2^4 #number of elements
n = 2^2 - 1 #polynomial order,

# for 64 total dof, K = 2^3, n = 2^3 - 1 is the break even point b/w sparse and full
# for K = 2^4, n = 2^2 - 1 sparse does better
# for K = 2^2, n = 2^4 - 1 full does better

println("The degrees of freedom are ")
println((n+1) * K)

# set domain parameters
L    = 2π
xmin = 0.0
xmax = L

# generate mesh variables
𝒢 = Mesh(K, n, xmin, xmax)

# generate internal variables
ι = Field1D(𝒢)

# set external parameters
ϰ = 1.0   #
α = 1.0   # parameter for solution, 1.0 is the example in the book
τ = 1.0  # penalty parameter
ε = external_params(ϰ, α)

# easy access
x  = 𝒢.x
u  = ι.u
u̇ = ι.u̇
q = copy(u)
dq = copy(ι.flux)

if periodic
    make_periodic1D!(𝒢.vmapP, ι.u)
end

f = 𝒢.M * sin.(α .* x) .* α^2
@. f *= 1 / 𝒢.rx
sol = -sin.(α * x)

params = (𝒢, ι, ε, periodic, q, dq, τ)

∇² = constructLaplacian(𝒢, periodic, τ)

∇² = Symmetric(∇²)
display(∇²)

#for plotting
theme(:juno)


s∇²  = sparse(∇²)
