include(pwd() * "/symbolics/dg_eval_rules.jl")
##
# Domain and Boundary
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)
# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);

K = 8      # Number of elements
n = 1      # Polynomial Order
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
α = 0.2; # Rusanov parameter
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
y_dg = Data(u0);
u̇ = Data(nothing);
u = Field(y_dg, field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
κ = 0.001 # Diffusivity Constant

# Burgers equation rhs
pde_equation = [
    Explicit(u̇, RK2) == Explicit(-∂xᴿ(u * u * 0.5), RK2)  + Implicit(κ * ∂xᶜ(∂xᶜ(u)), RK2),
]

pde_equation = [
    σ == ∂xᶜ(u),
    u̇ == -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(σ)
]

rhs =  -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u))
-∂xᴿ(u * u * 0.5) = first_order_terms(rhs)
κ * ∂xᶜ(σ) = aux_variables(rhs)

pde_meta_data = Dict("name" => "Burgers Equation", 
                    "method" => "discontinuous Galerkin")
pde_system = PDESystem(pde_equation,
                       Ω;
                       initial_condition=u0,
                       bcs=nothing,
                       metadata=pde_meta_data,)
##
# expr = :(u̇ = -∂xᴿ(u * u * 0.5)  + κ * ∂xᶜ(∂xᶜ(u));); 
# to change expr.args[1].args[2].args[2].args[2].args[1] = :∂xᶜ; eval(expr)
# ODE set up
Δx = mesh.x[2] - mesh.x[1]
dt = minimum([Δx^2 / κ * 0.05, abs(Δx / α)*0.05])
p = (pde_system, u, dt);

function dg_burgers!(v̇, v, p, t, α = true, β = false)
    # unpack parameters
    pde_system = p[1]
    u = p[2]
    u.data.data .= real.(v)
    v̇ .= α * compute(pde_system.equations[1].rhs) + β * v̇
    return v̇
end

rhs! = dg_burgers!
tspan = (0.0, 20.0)
# Define ODE problem
##
using TimeMachine
using DiffEqBase
prob = IncrementingODEProblem(rhs!, u0, tspan, p);
# Solve it
ode_method = LSRK144NiegemannDiehlBusch()
solve(prob, ode_method; dt=dt, adjustfinal=true);

# Plot it
##
theme(:juno)
plt = plot(
    x,
    u0,
    xlims=(Ω.a, Ω.b),
    ylims = (-1.1,1.1),
    marker = 3,
    leg = false,
)
display(plt)
