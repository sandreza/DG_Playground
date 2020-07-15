include(pwd()*"/symbolics" * "/dg_eval_rules.jl")

# Mesh Details
K = 20     # Number of elements
n = 1      # Polynomial Order
a = 0.0 # left endpoint of domain
b = 2π  # right endpoint of domain
mesh = Mesh(K, n, a, b, periodic = true) # Generate Uniform Periodic Mesh
x = mesh.x

# initial condition
u0 = @. exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
# DG Objects set up
α = 0.2; # Rusanov parameter
field_md = DGMetaData(mesh, nothing, nothing);
central = DGMetaData(mesh, u0, Rusanov(0.0));
rusanov = DGMetaData(mesh, u0, Rusanov(α));
y_dg = Data(u0);
u̇ = Data(nothing);
u = Field(y_dg, field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
κ = 0.001 # Diffusivity Constant

# Burgers equation rhs
pde_equation = [
    u̇ == -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u))
]
pde_meta_data = Dict("name" => "Burgers Equation", "method" => "discontinuous Galerkin")
pde_system = PDESystem(pde_equation,
                       mesh,
                       u0,
                       nothing,
                       pde_meta_data)

# expr = :(u̇ = -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u));); 
# to change expr.args[1].args[2].args[2].args[2].args[1] = :∂xᶜ; eval(expr)
# ODE set up
p = (u̇, u)

function dg_burgers!(v̇ , v, params, t)
    # unpack parameters
    u̇ = params[1]           
    u = params[2]
    u.data.data .= real.(v)
    v̇ .= eval(pde_system.equations[1].rhs)
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