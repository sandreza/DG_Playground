
include(pwd() * "/symbolics" * "/dg_eval_rules.jl")
include(pwd() * "/symbolics" * "/weno_functions.jl")
gr(size = (500,500))
for jjj in [1] ##[1, 2,3,4,5,6,7]
    for quadrature_rule in [true] # [true, false]
        for problemtype in ["burgers"] # ["burgers", "advection"]
n = jjj
problem = problemtype
inexact = quadrature_rule
information = "doing n="*string(n)*"_inexact="*string(inexact)*"_problem="*problem
println(information)

if inexact
    title = "diagonal mass"
else
    title = "block mass"
end

if problem == "burgers"
    T = 1.5 
elseif problem=="advection"
    c = 1.5
    T = 4π/c
else
    T = 1.5 
end

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
cfl = 0.3
K = 80   # Number of elements

r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
avg = Diagonal(zeros(n+1))
avg[1] = 1
cells = V * avg * inv(V)

mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
weights = sum(mesh.M, dims = 1)
Δx = x[2] - x[1]
if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end
b = Ω.b
a = Ω.a
ν = 1e-8
if problem == "advection"
    u0  = @. (tanh((x-3(b+a)/8)/ν)+1)*(tanh(-(x-5(b+a)/8)/ν)+1)/4 * (0.1*sin(40π/(b-a) * x) +2) # use initial condition for array
else
    u0  = @.sin(2π/(b-a) * x) + 0.5
end

α = 1*1.5
parameterlabel = "_α="*string(α)*"_"*title*"_p="*string(n)
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
u̇ = Data(u0);
u = Field(Data(u0), field_md);

∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
dt = minimum([abs(Δx / α) * cfl, abs(Δx / maximum(abs.(u0))) * cfl ])
dt = T / ceil(Int,T/dt)
numsteps = ceil(Int,T/dt)
if problem == "burgers"
    pde_equation = [
        u̇ == -∂xᴿ(u*u)*0.5,
    ]
    criteria = n*(mesh.x[2] - mesh.x[1])
elseif problem == "advection"
    pde_equation = [
        u̇ == -∂xᴿ(u*c)*0.5 ,
    ]
    criteria = n*100.0*(mesh.x[2] - mesh.x[1])
else
    pde_equation = []
end
equations = pde_equation

check = floor(Int,numsteps/20)
modby = check > 0 ? check : floor(Int,numsteps/4)
sol = []
push!(sol, copy(u.data.data))
for i in 1:numsteps   
    ssp3_step!(pde_equation, mesh, cells, dt, M = criteria)
    if (i%modby)==0
        push!(sol, copy(u.data.data))
    end
end
    if problem == "advection"
    ylims = (-0.2, 2.3)
else
    ylims = (-0.6, 1.6)
end

anim = @animate for i in eachindex(sol)

    v = sol[i]
    i1, v̅ =  troubled_cells(v, mesh, cells, M = criteria)
    plt = plot(mesh.x, v, 
    label = false, 
    linewidth =3,
    title= title * ", p = " * string(n),
    color = :blue, ylims = ylims)
    plot!(mesh.x[:,i1], v[:,i1], 
    label = false,
    linewidth =3, 
    color = :red, ylims = ylims)
end
gif(anim, pwd() * "/fig/weno_"*problem*"_animation"*parameterlabel*".gif")

fr = collect(range(-1,1,length=100))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)
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
    p1 = plot!(refine*x, refine*u.data.data,
     xlims = (3.8,4.0), 
     linewidth = 3,
     label = false)
    savefig(p1, pwd() * "/fig/weno_"*problem*"_ref_v_comp"*parameterlabel*"_weno.png")
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
    p1 = plot!(refine*x, 
    refine*u.data.data, 
    linewidth = 3,
    label = false)
    savefig(p1, pwd() * "/fig/weno_"*problem*"_ref_v_comp"*parameterlabel*".png")
end

end
end
end