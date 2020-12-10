include(pwd()*"/symbolics" * "/dg_eval_rules.jl")
include(pwd()*"/symbolics" * "/weno_functions.jl")
using JLD2
n = 4      # polynomial order
K = 80   # Number of elements
problem = "burgers" # "burgers" or "advection"
inexact = true # true or false quadrature_rule
weno = true
timescale = 1.0
information = "doing n="*string(n)*"_inexact="*string(inexact)*"_problem="*problem
println(information)
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
fr = collect(range(-1,1,length= 6*(n+1)))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)

if inexact
    title = "diagonal mass"
else
    title = "block mass"
end

if problem == "burgers"
    T = timescale * 1.5 
elseif problem=="advection"
    c = 1.5
    T = timescale * 4π/c
else
    T = 1.5 
end

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
cfl = 0.3


r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
avg = Diagonal(zeros(n+1))
avg[1] = 1
cells = V * avg * inv(V)

mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
weights = sum(mesh.M, dims = 1)
Δx = x[2] - x[1]
##
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

α = 1*1.5 / 5
parameterlabel = "_α="*string(α)* "_" * title * "_p=" * string(n)
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
u̇ = Data(u0);
u = Field(Data(u0), field_md);
ϕ = Field(Data(u0 .* u0), field_md);

∂xᴿ(a::AbstractExpression) = Gradient(a, rusanov);
dt = minimum([abs(Δx / α) * cfl, abs(Δx / maximum(abs.(u0))) * cfl ])
dt = T / ceil(Int,T/dt)
numsteps = ceil(Int,T/dt)
if problem == "burgers"
    pde_equation = [
        u̇ == -∂xᴿ(u*u*0.5),
    ]
    criteria = 0.01 * (Δx)^2
elseif problem == "advection"
    pde_equation = [
        u̇ == -∂xᴿ(c * u) ,
    ]
    criteria = 100 * (Δx)^2 
else
    pde_equation = []
end
equations = pde_equation

check = floor(Int,numsteps/60)
modby = check > 0 ? check : floor(Int,numsteps/20)
sol = []
push!(sol, copy(u.data.data))
for i in 1:numsteps   
    ssp3_step!(pde_equation, mesh, cells, dt, M = criteria, weno = weno)
    if (i%modby)==0
        push!(sol, copy(u.data.data))
    end
end
    if problem == "advection"
    ylims = (-0.2, 2.3)
else
    ylims = (-0.6, 1.6)
end
##
gr(size=(300,300))
theme(:juno)
for i in eachindex(sol)
    v = sol[i]
    i1, v̅ =  troubled_cells(v, mesh, cells, M = criteria)
        gr(size=(300,300))
    plt = plot(refine * mesh.x, refine * v, 
    label = false, 
    linewidth = 3,
    title= title * ", p = " * string(n),
    color = :blue, ylims = ylims)
    plot!(refine * mesh.x[:,i1], refine * v[:,i1], 
    label = false,
    linewidth =2, 
    color = :red, 
    ylims = ylims)
    sleep(0.1)
    display(plt)
end
##
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
##
gr(size = (300,300))
i = size(sol)[1]

xlims = (0, 2π)
xlims = (3.8, 4.0)
v = sol[i] 
i1, v̅ =  troubled_cells(v, mesh, cells, M = criteria)
plt = plot(refine * mesh.x, refine * v, 
label = false, 
linewidth = 3,
title= title * ", p = " * string(n),
color = :blue, xlims = xlims)
plot!(refine * mesh.x[:,i1], refine * v[:,i1], 
label = false,
linewidth = 3, 
color = :red,
)
sleep(0.1)
display(plt)
##
β = [smoothness_indicator(v[:,i], i, mesh) for i in eachindex(i1)]
cellindex = i1[argmax(β)]
cellindex = 52
vl, vr = extrapolate(v, mesh, cellindex)
plot(refine * mesh.x[:,cellindex], refine * v[:,cellindex], 
label = false,
linewidth =3, 
color = :red,
)
plot!(refine * mesh.x[:, cellindex .- 1 ], refine * v[:,cellindex .- 1], 
label = false,
linewidth =3, 
color = :blue,
)
plot!(refine * mesh.x[:,cellindex], refine * vl, 
label = false,
linewidth =3, 
color = :green,
)
plot!(refine * mesh.x[:, cellindex], refine * vr, 
label = false,
linewidth =3, 
color = :green,
)
plot!(refine * mesh.x[:,cellindex .+ 1], refine * v[:,cellindex .+ 1], 
label = false,
linewidth =3, 
color = :blue,
)
##
weights = sum(mesh.M, dims = 1) ./ 2
vl = vl .- weights * vl + v̅[:,cellindex]
vr = vr .- weights * vr + v̅[:,cellindex]
plot(refine * mesh.x[:,cellindex], refine * v[:,cellindex], 
label = false,
linewidth =3, 
color = :red,
)
plot!(refine * mesh.x[:, cellindex .- 1 ], refine * v[:,cellindex .- 1], 
label = false,
linewidth =3, 
color = :blue,
)
plot!(refine * mesh.x[:,cellindex], refine * vl, 
label = false,
linewidth =3, 
color = :green,
)
plot!(refine * mesh.x[:, cellindex], refine * vr, 
label = false,
linewidth =3, 
color = :green,
)
plot!(refine * mesh.x[:,cellindex .+ 1], refine * v[:,cellindex .+ 1], 
label = false,
linewidth =3, 
color = :blue,
)

##
β₀ =  smoothness_indicator(vl, cellindex, mesh)
β₁ =  smoothness_indicator(v[:, cellindex], cellindex, mesh)
β₂ =  smoothness_indicator(vr, cellindex, mesh)
(v[end, cellindex-1] -v[1, cellindex])^2 
(v[1, cellindex+1] -v[end, cellindex])^2 
vl, vr = extrapolate(v, mesh, cellindex)
vre = (vr[end]+v[end, cellindex])/2
vle = (vl[1]+v[1, cellindex])/2
vl = vl .- weights * vl + v̅[:,cellindex]
vr = vr .- weights * vr + v̅[:,cellindex]

β = [β₀ β₁ β₂]
ϵ = 1e-6
allw = 1 ./ ((ϵ .+ β) .^ 2)
γ1 = 0.998
γ0 = (1-γ1)/2
γ2 = γ0
sumw = γ0*allw[1] + γ1*allw[2] + γ2*allw[3]
w0 = γ0*allw[1] ./ sumw
w1 = γ1*allw[2] ./ sumw
w2 = γ2*allw[3] ./ sumw

reconstruction = w0*vl + w1*v[:, cellindex] + w2*vr
smoothness_indicator(reconstruction, cellindex, mesh)
plot!(refine * mesh.x[:,cellindex], refine * reconstruction, 
label = false,
linewidth =3, 
color = :purple,
)
#=
## plot
## Step 1 Identify troubled cells
gr(size=(1000,1000))
cellindex = i1[end-2]
cellindexr = adjust(cellindex+1, mesh)
cellindexl = adjust(cellindex-1, mesh)
vr, vl = extrapolate(v, mesh, cellindex)
fr = collect(range(-1,1,length=100))
fV = vandermonde(fr, 0, 0, n)
refine = fV * inv(V)

identifyp = plot(refine*x[:, cellindex], 
refine*v[:,cellindex],
color = :red, label = "troubled cell",
legend = :right,
linewidth = 3)
plot!(refine*x[:, cellindexl],
 refine*v[:, cellindexl],
 color = :blue,
  label = "neighborl",
  linewidth = 3)
plot!(refine*x[:, cellindexr],
 refine*v[:, cellindexr],
color = :blue,
label = "neighborr",
title = "identify",
linewidth = 3)


extrapolatep = plot(refine*x[:, cellindex], 
refine*v[:,cellindex],
color = :red, label = "troubled cell",
legend = :right,
linewidth = 3)
plot!(refine*x[:, cellindexl],
 refine*v[:, cellindexl],
 color = :blue,
  label = "neighborl",
  linewidth = 3)
plot!(refine*x[:, cellindexr],
 refine*v[:, cellindexr],
color = :blue,
label = "neighborr",
linewidth = 3)
plot!(refine*x[:, cellindex], 
refine*vr, 
color = :green,
label = "extrapolate r",
linewidth = 3)
plot!(refine*x[:, cellindex], 
refine*vl, 
color = :green,
label = "extrapolate l",
title = "extrapolate",
linewidth = 3)
# Step 2 Adjust
weights = sum(mesh.M, dims = 1) ./ 2
v̅  = weights * v[:, cellindex]
vr = vr .- weights*vr .+ v̅ 
vl = vl .- weights*vl .+ v̅ 
adjustp = plot(refine*x[:, cellindex], 
refine*v[:,cellindex],
color = :red, label = "troubled cell",
legend = :topright,
linewidth = 3)
plot!(refine*x[:, cellindexl],
 refine*v[:, cellindexl],
 color = :blue,
  label = "neighborl",
  linewidth = 3)
plot!(refine*x[:, cellindexr],
 refine*v[:, cellindexr],
color = :blue,
label = "neighborr",
linewidth = 3)
plot!(refine*x[:, cellindex], 
refine*vr, 
color = :green,
label = "extrapolate r",
linewidth = 3)
plot!(refine*x[:, cellindex], 
refine*vl, 
color = :green,
label = "extrapolate l",
title = "adjust",
linewidth = 3)
# Step 3: Calculate smoothness
β₀ =  smoothness_indicator(vl, cellindex, mesh)
β₁ =  smoothness_indicator(v[:, cellindex], cellindex, mesh)
β₂ =  smoothness_indicator(vr, cellindex, mesh)
β = [β₀ β₁ β₂]
ϵ = 1e-6
allw = 1 ./ ((ϵ .+ β) .^ 2)
γ1 = 0.95 #998
γ0 = (1-γ1)/2
γ2 = γ0
sumw = γ0*allw[1] + γ1*allw[2] + γ2*allw[3]
w0 = γ0*allw[1] ./ sumw
w1 = γ1*allw[2] ./ sumw
w2 = γ2*allw[3] ./ sumw
reconstruct = w0*vl + w1*v[:, cellindex] + w2*vr

reconstructp = plot(refine*x[:, cellindex], 
refine*v[:,cellindex],
color = :red, label = "troubled cell",
legend = :topright,
linewidth = 3)
plot!(refine*x[:, cellindex],
refine*reconstruct,
color = :green,
label = "reconstructed",
title = "reconstruct",
linewidth = 3)
plot!(refine*x[:, cellindexl],
 refine*v[:, cellindexl],
 color = :blue,
  label = "neighborl",
  linewidth = 3)
plot!(refine*x[:, cellindexr],
 refine*v[:, cellindexr],
color = :blue,
label = "neighborr",
linewidth = 3)

gr(size=(1000,1000))
pall = plot(identifyp, extrapolatep, adjustp, reconstructp)
savefig(pall, pwd() * "/fig/"*"weno_procedure.png")
##
=#
##
gr(size = (500,500))
expand = 0.11
shock = (3.8-expand,4.0+expand)
method = u.metadata.method
data = u.data.data
dlims = (minimum((∫dV + ∫dA)), maximum((∫dV + ∫dA)))
∫dV = compute_volume_terms(data, mesh)
∫dA = compute_surface_terms(mesh, data, data, Rusanov(0.0))
p1 = plot(refine * x, refine * u.data.data, xlims = shock, label = false, title = "function"  * " p="*string(n))
p2 = plot(refine * x, refine * ∫dV, xlims = shock, label= false, title = "\"volume\" derivative", ylims = dlims)
p3 = plot(refine * x, refine * ∫dA, xlims = shock, label= false, title = "\"surface\" derivative", ylims = dlims)
p4 = plot(refine * x, refine * (∫dV + ∫dA), xlims = shock, label = false, title = "total derivative", ylims = dlims)
plot(p1,p2,p3,p4)
