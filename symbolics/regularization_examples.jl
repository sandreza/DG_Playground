include(pwd()*"/symbolics" * "/dg_eval_rules.jl")
using DifferentialEquations

    println(n, problem, inexact)
# Initial Condition
if problem == "burgers"
    u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
    T = 1.5 
elseif problem=="advection"
    u⁰(x, a, b; ν = 0.0001) = (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4 * (0.1*sin(40π/(b-a) * x) +2)
    c = 1.5
    T = 4π/c
else
    u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
    T = 1.5 
end

if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end

u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
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
##
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
##
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5, Feagin14()
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it
##
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
##
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
display(plot(plt,plt2))
savefig(plot(plt,plt2), pwd() * "/fig/"*problem*"_last_step"*parameterlabel*".png")
