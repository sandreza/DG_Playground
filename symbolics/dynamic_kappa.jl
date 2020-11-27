include(pwd()*"/symbolics" * "/dg_eval_rules.jl")
# Domain and Boundary
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
#u⁰(x, a, b; ν = 0.0001) = (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4 * (0.1*sin(40π/(b-a) * x) +2)
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
c = 1.5
T = 1.5 # 4pi/1.5# advection = 4π/1.5, 1.5 for burgers
inexact = true
cfl = 0.3
K = 80   # Number of elements
n = 3    # Polynomial Order
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)

flux_limiter = false
exact_nonlinear = false
timestepfilter = false
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
# mesh.D .= mesh.D * linearfilter
Δx = x[2] - x[1]
plot(x, u⁰.(x, Ω.a,  Ω.b))
if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end

##

u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
α = 1.5/2#1.5/2; # Rusanov parameter
Δu = 0.5 / n# regularization parameter
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
u̇ = Data(u0);
u = Field(Data(u0), field_md);
u² = Field(Data(u0 .* u0), field_md);
κbase = minimum([1.5e-8, Δx * α ]) # Diffusivity Constant
κ0 = x .* 0 .+ κbase
kmax = Δx * α 
κ̇ = Data(κ0);
κ = Field(Data(κ0) , field_md);
∂xᶜ(a::AbstractExpression) = Gradient(a, central);
∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
dt = minimum([Δx^2 / maximum(kmax) * cfl, abs(Δx / α) * cfl, abs(Δx / maximum(abs.(u0))) * cfl ])
λ = 0.5/(dt) # 1.0 / dt makes it change immediately every timstep

# minimum grid spacing  / 2*expected linfinity
smoothness = tanh((∂xᶜ(u)*(Δx/(Δu+eps(1.0))))^2) # 0 for smooth, 1 for unsmooth
# Burgers equation rhs, 
pde_equation = [
    κ̇ == -λ*(κ + (-κbase + 1 * (-kmax+κbase) * smoothness)) ,
    u̇ == -∂xᴿ(u²)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
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
    pde_system.equations[1].lhs.data .=  real.(v[1:(n+1),:])
    pde_system.equations[2].lhs.data .= real.(v[(n+2):end,:])
    l_u = pde_system.equations[2].lhs.data
    u².data.data .= l_u .* l_u
    # modify for exact nonlinear
    if (n==1) & exact_nonlinear
        v̇[1:(n+1),:] .= compute(pde_system.equations[1].rhs)
        ω1 = 5/6
        ω2 = 1/3  
        if inexact
            ω1 = 1/2
            ω2 = 1/3
        end
        ω3 = 1 - ω1 - ω2
        for i in 1:K
            u².data.data[1,i]  =  ω1*l_u[1,i]*l_u[1,i] + ω2*l_u[1,i]*l_u[end,i] + ω3*l_u[end,i]*l_u[end,i]
            u².data.data[end,i] = ω3*l_u[1,i]*l_u[1,i] + ω2*l_u[1,i]*l_u[end,i] + ω1*l_u[end,i]*l_u[end,i]
        end
        # the above modifies the volume nodes, the surface nodes stay with the same estimate
        ∫dA1 = -0.5 * compute_surface_terms(u².metadata.mesh, l_u .* l_u, l_u, Rusanov(α))
        ∫dA2 = -0.5 * compute_surface_terms(u².metadata.mesh, u².data.data, l_u, Rusanov(α))
        v̇[(n+2):end,:] .= compute(pde_system.equations[2].rhs) + ∫dA1 - ∫dA2
    else
            if timestepfilter
        filteredthing = u².data.data
        spectrum = inv(V)*filteredthing
        if n > 2
            for j in 1:K
                for jj in 3:(n+1)
                    criteria = (abs.(spectrum[jj-1,j]) + abs.(spectrum[jj-2,j]))/2
                    if abs.(spectrum[jj,j]) > criteria
                        spectrum[jj,j] *= (0.9)^jj
                    end
                end
            end
        end
        u².data.data .= V * spectrum
    end
        v̇[1:(n+1),:] .= compute(pde_system.equations[1].rhs)
        v̇[(n+2):end,:] .= compute(pde_system.equations[2].rhs)
        if flux_limiter
            diagfilter = ones(n+1)
            diagfilter[1] = 1
            diagfilter[2] = 1.0
            if n+1>=3
                diagfilter[3:end] .= 1.0
            end
            filter = Diagonal(diagfilter)
            filter = V * filter * inv(V)
            ruf =  filter * (l_u .* l_u)
            ∫dA1 = -0.5 * compute_surface_terms(u².metadata.mesh, ruf, l_u, Rusanov(α))
            ∫dA2 = -0.5 * compute_surface_terms(u².metadata.mesh, u².data.data, l_u, Rusanov(α))
            v̇[(n+2):end,:] .+= ∫dA1 - ∫dA2
        end
    end

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
ode_method = Heun() # Heun(), RK4, Tsit5, Feagin14()
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
gr(size = (300,300))
#filtersolution
mesh = create_mesh(Ω, elements = K, polynomial_order =  n)
diagfilter = ones(n+1)
diagfilter[1] = 1
if n+1>=3
    diagfilter[3:end] = exp.(-0.01*(collect(3:1:n+1) .- 2.1 ))
end
filter = Diagonal(diagfilter)
linearfilter = V * filter * inv(V)
fr = collect(range(-1,1,length=100))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)
# anim = @animate  
if inexact
    title = "Inexact"
else
    title = "Exact"
end
for i in indices
    ylims = (minimum(sol.u[1])-0.1*maximum(sol.u[1]), maximum(sol.u[1]) + 0.1*maximum(sol.u[1]))
    plt = plot(refine*x, 
    refine*linearfilter * real.(sol.u[i])[(n+2):end,:], 
    xlims=(Ω.a, Ω.b), ylims = ylims,  
    linewidth = 2.0, 
    leg = false,
    title= title)
    plot!(x,  real.(sol.u[1])[(n+2):end,:], xlims = (Ω.a, Ω.b), ylims = ylims, color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    #plt = plot(x, real.(sol.u[i])[1:(n+1),:],ylims = (0, kmax), xlims = (Ω.a, Ω.b), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    data = sol.u[i][(n+2):end,:]
    # Legendre Mode Amplitudes
    p = []
    push!(p,scatter(log10.(abs.(inv(V) * data[:,1])), label = false ))
    # label = "element " * string(1)) 
    for j in 2:K
        push!(p,scatter(log10.(abs.(inv(V) * data[:,j])), label = false ))
    end
    plt2 = plot(p[1:K]...)
    display(plt)
    #display(plt2)
    sleep(0.05)
end
##
expand = 0.0
shock = (3.8-expand,4.0+expand)
method = u.metadata.method
data = u.data.data
∫dV = compute_volume_terms(data, mesh)
∫dA = compute_surface_terms(mesh, data, data, Rusanov(0.0))
p1 = plot(refine * x, refine * u.data.data, xlims = shock, label = false, title = title  * " p="*string(n))
p2 = plot(refine * x, refine * ∫dV, xlims = shock, label= false, title = "volume derivative")
p3 = plot(refine * x, refine * ∫dA, xlims = shock, label= false, title = "surface derivative")
p4 = plot(refine * x, refine * (∫dV + ∫dA), xlims = shock, label = false, title = "total derivative")
plot(p1,p2,p3,p4)
##
#=
tmp = maximum(real.(sol.u[end])[1:(n+1),:])
plot(x, real.(sol.u[end])[1:(n+1),:], xlims = (Ω.a, Ω.b), ylims = (0.00, tmp), leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
##

if inexact
    gif(anim, "burgers_inexact_2.gif")
else
    gif(anim, "burgers_inexact_1.gif")
end


tmp = real.(sol.u[end])[(n+2):end,:]

tmp[1,:] - tmp[end,:]
=#