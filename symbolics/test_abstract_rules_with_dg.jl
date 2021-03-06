include(pwd()*"/symbolics" * "/dg_eval_rules.jl")

# Domain and Boundary
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x, a, b; ν = 0.1) = (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
c = 2π/10
inexact = true
cfl = 0.25
K = 80     # Number of elements
n = 4    # Polynomial Order
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
α = 10.0; # Rusanov parameter
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
y_dg = Data(u0);
u̇ = Data(nothing);
u = Field(y_dg, field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
κ = 1e-1 # Diffusivity Constant

# Burgers equation rhs
pde_equation = [
    u̇ == -∂xᴿ(u * u)*0.5  + κ * ∂xᶜ(∂xᶜ(u)),
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
p = (pde_system, u);

function dg_burgers!(v̇ , v, params, t)
    # unpack parameters
    pde_system = params[1]
    u = params[2]
    u.data.data .= real.(v)
    v̇ .= compute(pde_system.equations[1].rhs)
    return nothing
end

rhs! = dg_burgers!
tspan = (0.0, 1.5)

# Define ODE problem
ode_problem = (rhs!, u0, tspan, p);
##
using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = RK4() # Heun(), RK4, Tsit5
Δx = mesh.x[2] - mesh.x[1]
dt = minimum([Δx^2 / κ * cfl, abs(Δx / α) * cfl, abs(Δx / maximum(abs.(u0))) * cfl ])
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
    plt = plot(x, real.(sol.u[i]), color = :blue, xlims=(Ω.a, Ω.b), ylims = ylims, marker = 3,  leg = false)
    plot!(x, real.(sol.u[1]), xlims = (Ω.a, Ω.b), ylims = ylims, color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    sleep(0.05)
end
reference = sol.u[end]
if inexact
    gif(anim, "burgers_inexact_2.gif")
else
    gif(anim, "burgers_inexact_1.gif")
end
##
plot(x, real.(sol.u[end]), xlims=(Ω.a, Ω.b), ylims = (-1.1,1.1), marker = 3,  leg = false)
plot!(ref_grid, ref_sol, xlims = (a, b), ylims = (-1.1,1.1), color = "blue", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, line = 3, label = "Reference Solution")