# n is the Polynomial Order
using JLD2, DifferentialEquations
include(pwd()*"/symbolics" * "/dg_eval_rules.jl")

for jjj in [1, 2,3,4,5,6,7]
    for quadrature_rule in [true, false]
        for problemtype in ["burgers", "advection"]
n = jjj
problem = problemtype
inexact = quadrature_rule

information = "doing n="*string(n)*"_inexact="*string(inexact)*"_problem="*problem
println(information)

# all in the loop   
cfl = 0.3
K = 80   # Number of elements


Ω  = IntervalDomain(0, 2π, periodic = true)
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
x = mesh.x
Δx = mesh.x[2] - mesh.x[1]

α = 1.5/2    # Rusanov Wavespeed
Δu = 1.5 / n # regularization parameter
kmax = Δx * α 
κbase = 1e-8 # Diffusivity Constant

dt = minimum([Δx^2 / maximum(kmax) * cfl, abs(Δx / α) * cfl, abs(Δx / 2.0) * cfl ])
if inexact
    title = "diagonal mass"
else
    title = "block mass"
end

λ = 0.5/(dt) # 1.0 / dt makes it change immediately every timstep
parameterlabel = "_α="*string(α)*"_"*title*"_p="*string(n)


println(n, problem, inexact)
# Initial Condition


if problem == "burgers"
    T = 1.5 
elseif problem=="advection"
    c = 1.5
    T = 4π/c
else
    T = 1.5 
end

if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end
ν = 0.001
b = Ω.b
a = Ω.a

if problem == "advection"
    u0  = @. (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4 * (0.1*sin(40π/(b-a) * x) +2) # use initial condition for array
else
    u0  = @.sin(2π/(b-a) * x) + 0.5
end

field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata

u̇ = Data(u0);
u = Field(Data(u0), field_md);
κ0 = x .* 0 .+ κbase
κ̇ = Data(κ0);
κ = Field(Data(κ0) , field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);

# minimum grid spacing  / 2*expected linfinity
smoothness = tanh((∂xᶜ(u)*(Δx/(Δu+eps(1.0))))^2) # 0 for smooth, 1 for unsmooth
# Burgers equation rhs, 
if problem == "burgers"
    pde_equation = [
        κ̇ == -λ*(κ + (-κbase + 1 * (-kmax+κbase) * smoothness)) ,
        u̇ == -∂xᴿ(u*u)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
    ]
elseif problem == "advection"
    pde_equation = [
        κ̇ == -λ*(κ + (-κbase + 1 * (-kmax+κbase) * smoothness)) ,
        u̇ == -∂xᴿ(u*c)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
    ]
else
    pde_equation = [
        κ̇ == -λ*(κ + (-κbase + 1 * (-kmax+κbase) * smoothness)) ,
        u̇ == -∂xᴿ(u*u)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
    ]
end

pde_meta_data = Dict("name" => "Burgers Equation", "method" => "discontinuous Galerkin")
pde_system = PDESystem(pde_equation,
                       Ω;
                       initial_condition=u0,
                       bcs=nothing,
                       metadata=pde_meta_data)
p = (pde_system, u, κ, n);
function dg_burgers!(v̇ , v, params, t)
    # unpack parameters
    pde_system = params[1]
    u = params[2]
    κ = params[3]
    n = params[4]
    pde_system.equations[1].lhs.data .=  real.(v[1:(n+1),:])
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

prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5, Feagin14()
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it

theme(:juno)
nt =  length(sol.t)
num = 40 # Number of Frames
stp = floor(Int, nt/num)
num = floor(Int, nt/stp)
indices = stp * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)
gr(size = (500,500))

fr = collect(range(-1,1,length=100))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)
# anim = @animate  

anim = @animate for i in indices
    ylims = (minimum(sol.u[1])-0.1*maximum(sol.u[1]), maximum(sol.u[1]) + 0.1*maximum(sol.u[1]))
    plt = plot(refine*x, 
    refine*real.(sol.u[i])[(n+2):end,:], 
    xlims=(Ω.a, Ω.b), ylims = ylims,  
    linewidth = 1.0, 
    leg = false,
    title= title * ", p = " * string(n),
    xlabel = "x",
    ylabel = "u")
    if problem=="advection"
            plt = plot!(refine*x, 
        refine*real.(sol.u[1])[(n+2):end,:], 
        xlims=(Ω.a, Ω.b), ylims = ylims,  
        linewidth = 1.0, 
        leg = false,
        title= title * ", p = " * string(n),
        xlabel = "x",
        ylabel = "u",
        color = :red)
    end

    plt2 = plot(x, log10.(real.(sol.u[i])[1:(n+1),:]),
    ylims = (log10(κbase)-0.1, log10(kmax)), xlims = (Ω.a, Ω.b), 
    color = "red", 
    leg = false, grid = true, 
    gridstyle = :dash, gridalpha = 0.25, 
    framestyle = :box,
    xlabel = "x",
    ylabel = "log10(κ)",
    title = "smoothness indicator")
    plot(plt,plt2)
end
gif(anim, pwd() * "/fig/"*problem*"_animation"*parameterlabel*".gif")

i = nt
ylims = (minimum(sol.u[1])-0.1*maximum(sol.u[1]), maximum(sol.u[1]) + 0.1*maximum(sol.u[1]))
plt = plot(refine*x, 
refine*real.(sol.u[i])[(n+2):end,:], 
xlims=(Ω.a, Ω.b), ylims = ylims,  
linewidth = 1.0, 
leg = false,
title= title * ", p = " * string(n),
xlabel = "x",
ylabel = "u")

plt2 = plot(x, log10.(real.(sol.u[i])[1:(n+1),:]),
ylims = (log10(κbase)-0.1, log10(kmax)), xlims = (Ω.a, Ω.b), 
color = "red", 
leg = false, grid = true, 
gridstyle = :dash, gridalpha = 0.25, 
framestyle = :box,
xlabel = "x",
ylabel = "log10(κ)",
title = "smoothness indicator")

savefig(plot(plt,plt2), pwd() * "/fig/"*problem*"_last_step"*parameterlabel*".png")

if problem == "burgers"
    refsolution = jldopen("reference.jld2")
    refx = refsolution["x"]
    refu = refsolution["u.data.data"]
    plot(refx, refu, 
    linewidth = 3, 
    xlims = (3.8, 4.0), 
    label = false, 
    color = :blue,
    title = "reference vs computed, "*title*", p="*string(n),
    ylabel = "u",
    xlabel = "x")
    p1 = plot!(refine*x, refine*u.data.data, xlims = (3.8,4.0), label = false)
    savefig(p1, pwd() * "/fig/"*problem*"_ref_v_comp"*parameterlabel*".png")
else
    refx = refine * x
    refu  = @. (tanh((refx-3(b+a)/8)/ν)+1)*(tanh(-(refx-5(b+a)/8)/ν)+1)/4 * (0.1*sin(40π/(b-a) * refx) +2)
    plot(refine * x, refu, 
    linewidth = 3, 
    label = false, 
    color = :blue,
    title = "reference vs computed, "*title*", p="*string(n),
    ylabel = "u",
    xlabel = "x")
    p1 = plot!(refine*x, refine*u.data.data, label = false)
    savefig(p1, pwd() * "/fig/"*problem*"_ref_v_comp"*parameterlabel*".png")
end

        end
    end
end
