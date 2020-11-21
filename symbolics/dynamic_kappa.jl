include(pwd()*"/symbolics" * "/dg_eval_rules.jl")

# Domain and Boundary
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x, a, b; ν = 0.001) = (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
c = 2π/10
T = 1.5 #20 # 1.5 for burgers
inexact = true
cfl = 0.25
K = 100     # Number of elements
n = 1    # Polynomial Order
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
plot(x, u⁰.(x, Ω.a,  Ω.b))
if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end

u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
α = 0.1; # Rusanov parameter
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
u̇ = Data(u0);
u = Field(Data(u0), field_md);
κbase = 1e-4 # Diffusivity Constant
κ0 = x .* 0 .+ κbase
kmax = 1e-1
κ̇ = Data(κ0);
κ = Field(Data(κ0) , field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
dt = minimum([Δx^2 / maximum(kmax) * cfl, abs(Δx / α) * cfl, abs(Δx / maximum(abs.(u0))) * cfl ])
λ = 0.1/(dt) # 1.0 / dt makes it change immediately every timstep

smoothness = tanh(∂xᶜ(u)^2*0.01) # 0 for smooth, 1 for unsmooth
# Burgers equation rhs, 
pde_equation = [
    κ̇ == -λ*(κ + (-κbase + 1 * (-kmax) * smoothness)),
    u̇ == -∂xᴿ(u * u)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
]

pde_meta_data = Dict("name" => "Burgers Equation", "method" => "discontinuous Galerkin")
pde_system = PDESystem(pde_equation,
                       Ω;
                       initial_condition=u0,
                       bcs=nothing,
                       metadata=pde_meta_data)

##
# expr = :(u̇ = -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u));); 
# to change expr.args[1].args[2].args[2].args[2].args[1] = :∂xᶜ; eval(expr)
# ODE set up
p = (pde_system, u, κ, n);

function dg_burgers!(v̇ , v, params, t)
    # unpack parameters
    pde_system = params[1]
    u = params[2]
    κ = params[3]
    n = params[4]
    pde_system.equations[1].lhs.data .= real.(v[1:(n+1),:])
    pde_system.equations[2].lhs.data .= real.(v[(n+2):end,:])
    v̇[1:(n+1),:] .= compute(pde_system.equations[1].rhs)
    v̇[(n+2):end,:] .= compute(pde_system.equations[2].rhs)
    return nothing
end

rhs! = dg_burgers!
tspan = (0.0, T)

# Define ODE problem
v = vcat(copy(κ0), copy(u0))
v̇ = (copy(κ0), copy(u0))
ode_problem = (rhs!, v, tspan, p);
##
using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = RK4() # Heun(), RK4, Tsit5
Δx = mesh.x[2] - mesh.x[1]
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it
##
theme(:juno)
nt =  length(sol.t)
num = 10 # Number of Frames
stp = floor(Int, nt/num)
num = floor(Int, nt/stp)
indices = stp * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
anim = @animate  for i in indices
    ylims = (minimum(sol.u[1])-0.1*maximum(sol.u[1]), maximum(sol.u[1]) + 0.1*maximum(sol.u[1]))
    plt = plot(x, real.(sol.u[i])[(n+2):end,:], color = :blue, xlims=(Ω.a, Ω.b), ylims = ylims, marker = 3,  leg = false)
    plot!(x, real.(sol.u[1])[(n+2):end,:], xlims = (Ω.a, Ω.b), ylims = ylims, color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    sleep(0.05)
end
##
tmp = maximum(real.(sol.u[end])[1:(n+1),:])
plot(x, real.(sol.u[end])[1:(n+1),:], xlims = (Ω.a, Ω.b), ylims = (0.00, tmp), leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
##
if inexact
    gif(anim, "burgers_inexact_2.gif")
else
    gif(anim, "burgers_inexact_1.gif")
end

##

tmp = real.(sol.u[end])[(n+2):end,:]

tmp[1,:] - tmp[end,:]