using SparseArrays, BenchmarkTools, Plots
include(pwd() * "/symbolics/abstract_core.jl")
include(pwd() * "/src" * "/HesthavenWarburton" * "/utils.jl")
include(pwd() * "/src" * "/HesthavenWarburton" * "/mesh.jl")
##

function compute_volume_terms(data::AbstractArray, mesh::Mesh)
    q = mesh.D * data
    @. q *= mesh.rx
    return q
end

struct Rusanov{𝒯}
    α::𝒯
end
##
function compute_surface_terms(mesh::AbstractMesh, data, state::AbstractArray, method::Rusanov{𝒯}) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (data[mesh.vmapM] + data[mesh.vmapP]), (mesh.nFP * mesh.nFaces, mesh.K ))

    # Include factor of 2 for the weak-strong form
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusanov
    @. diffs[:] += method.α * mesh.normals[:] .* (state[mesh.vmapM] - state[mesh.vmapP]) / 2.0
    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= data[mesh.vmapM]
    # Compute Lift Operator
    lifted =  mesh.lift * (mesh.fscale .* mesh.normals .* diffs)
    return lifted
end

struct Gradient{𝒯, 𝒰} <: AbstractExpression
    operand::𝒯
    metadata::𝒰
end

function dg_derivative(mesh, data, state, method)
    ∫dV = compute_volume_terms(data, mesh)
    ∫dA = compute_surface_terms(mesh, data, state, method)
    return ∫dV .+ ∫dA
end

struct DGMetaData{𝒮, 𝒯, 𝒰} 
    mesh::𝒮
    state::𝒯
    method::𝒰
end
##
dg_derivative(y::AbstractArray, md) = dg_derivative(md.mesh, y, md.state, md.method)
dg_derivative(y::AbstractData, md) = dg_derivative(md.mesh, y.data, md.state, md.method)
function eval(e::Gradient{𝒯, 𝒰}) where {𝒯, 𝒰 <: DGMetaData}
    return dg_derivative(eval(e.operand), e.metadata)
end

##
K = 20     # Number of elements
n = 1      # Polynomial Order
a = 0.0 # left endpoint of domain
b = 2π  # right endpoint of domain
mesh = Mesh(K, n, a, b, periodic = true) # Generate Uniform Periodic Mesh
x = mesh.x
D = mesh.D
volume = mesh.rx


# initial condition
u0 = @. exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
α = 0.2;
field_md = DGMetaData(mesh, nothing, nothing);
central = DGMetaData(mesh, u0, Rusanov(0.0));
rusanov = DGMetaData(mesh, u0, Rusanov(α));
y_dg = Data(u0);
u = Field(y_dg, field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
κ = 0.001
# Burgers equation rhs
u̇ = -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u));
p = (u̇, u)

function dg_burgers!(v̇ , v, params, t)
    # unpack params
    u̇ = params[1]           
    u = params[2]
    u.data.data .= real.(v)
    v̇ .= eval(u̇)
    return nothing
end

rhs! = dg_burgers!
tspan = (0.0, 20.0)

# Define ODE problem
ode_problem = (rhs!, u0, tspan, p);

##
using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5
Δx = mesh.x[2] - mesh.x[1]
dt = minimum([Δx^2 / κ * 0.05, abs(Δx / α)*0.05]) 
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it
##
theme(:juno)
nt = length(sol.t)
num = 40 # Number of Frames
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
for i in indices
    plt = plot(x, real.(sol.u[i]), xlims=(a, b), ylims = (-1.1,1.1), marker = 3, color = "green",   leg = false)
    plot!(x, real.(sol.u[1]), xlims = (a, b), ylims = (-1.1,1.1), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    sleep(0.1)
end