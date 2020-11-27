include(pwd()*"/symbolics" * "/dg_eval_rules.jl")
# Domain and Boundary
Ω  = IntervalDomain(0, 2π, periodic = true)
∂Ω = ∂(Ω)

# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
T = 1.5
inexact = true
cfl = 0.3
K = 80   # Number of elements
n = 4    # Polynomial Order
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
flux_limiter = false
exact_nonlinear = false
timestepfilter = false
mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
weights = sum(mesh.M, dims = 1)
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
κbase = minimum([1.5e-8, Δx * α ]) # Diffusivity Constant
κ0 = x .* 0 .+ κbase
kmax = Δx * α * 1 /10 + κbase
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
    u̇ == -∂xᴿ(u*u)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
]

p = (pde_equation, u, κ, n);

function dg_burgers!(v̇ , v, params, t)
    # unpack parameters
    equations = params[1]
    u = params[2]
    κ = params[3]
    n = params[4]
    equations[1].lhs.data .=  real.(v[1:(n+1),:])
    equations[2].lhs.data .=  real.(v[(n+2):end,:])
    v̇[1:(n+1),:]   .= compute(equations[1].rhs)
    v̇[(n+2):end,:] .= compute(equations[2].rhs)
    return nothing
end

rhs! = dg_burgers!
tspan = (0.0, T)

# Define ODE problem
v = vcat(copy(κ0), copy(u0))
v̇ = (copy(κ0), copy(u0))
ode_problem = (rhs!, v, tspan, p);

using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5, Feagin14()
sol  = solve(prob, ode_method, dt=dt, adaptive = false);

# Plot it

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
V = vandermonde(r, 0, 0, n)
avg = Diagonal(zeros(n+1))
avg[1] = 1
cells = V * avg * inv(V)
v = u.data.data
v̅ = cells * v
ṽ = v̅[1,:] -  v[1,:]
ṽ̃ = v[end, :] - v̅[end,:]
vtb = reshape(v̅[mesh.vmapP] - v̅[mesh.vmapM], (2,80))
Δ₊v = vtb[end,:]
Δ₋v = -vtb[1,:]
function minmod(a)
    if prod(sign.(a) .== sign(a[1]))
        return sign.(a[1])*minimum(abs.(a))
    else
        return -0
    end
end
# checks two things:
# 1
# If the jump in the cell average between neighbors
# is smaller than the jump between a DG cell average and
# its endpoints, then the cell is troubled
# 2
# If there is a change in sign between the jumps 
# the cell is trouble
mminmod(a; h=0.1, M=0.1) = abs.(a[1]) < h*M^2 ? abs.(a[1]) : minmod(a) 
troubled_i = collect(1:K)
minmods = [minmod([ṽ[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
troubled = minmods .!= ṽ
minmods = [minmod([ṽ̃[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
troubled2 = minmods .!= ṽ̃
makeitdouble = troubled .| troubled2
i1 = troubled_i[makeitdouble]
i2 = troubled_i[broadcast(~, makeitdouble)]
plt = plot(x[:,i1], v[:,i1], label = false, color = :red)
plot!(x[:,i2], v[:,i2], label = false, color = :blue)

for i in indices
    v = real.(sol.u[i])[(n+2):end,:]
    v̅ = cells * v
    ṽ = v̅[1,:] -  v[1,:]
    ṽ̃ = v[end, :] - v̅[end,:]
    vtb = reshape(v̅[mesh.vmapP] - v̅[mesh.vmapM], (2,80))
    Δ₊v = vtb[end,:]
    Δ₋v = -vtb[1,:]
    minmods = [minmod([ṽ[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
    troubled = minmods .!= ṽ
    minmods = [minmod([ṽ̃[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
    troubled2 = minmods .!= ṽ̃
    makeitdouble = troubled .| troubled2
    i1 = troubled_i[makeitdouble]
    i2 = troubled_i[broadcast(~, makeitdouble)]
    plt = plot(x[:,i1], v[:,i1], label = false, color = :red, ylims = (-0.6, 1.6))
    plot!(x[:,i2], v[:,i2], label = false, color = :blue)
    display(plt)
    sleep(0.05)
end

##
shock = (3.5, 4.4)
kmax = Δx * 0.75
p1 = plot(x, v .* v .* 0.5 - compute(kmax * ∂xᶜ(u)) , label = false, ylims = (-0.0, 2.25) , xlims = shock)
p2 = plot(x,v .* v .* 0.5 , label = false, ylims = (-0.0, 2.26) , xlims = shock)
plot(p1,p2)
##


v̅ = cells * v
ṽ = v-v̅

w1 = 0.3
w2 = 0.4
w3 = 1-w1-w2
plot(x[:,i1], ṽ[:,i1], label = false, color = :red)
plot!(x[:,i1], w1 .* ṽ[:,i1 .- 1] + w2 .* ṽ[:,i1]+ w3 .* ṽ[:,i1 .+ 1], label = false, color = :green)
#plot!(x[:,i2], ṽ[:,i2], label = false, color = :blue)
##
# smoothness
plot(x, mesh.D *(ṽ))
smoothness = zeros(K)
Δxⱼ = reshape(mesh.x[end,:]-mesh.x[1,:], (1,K))
deriv = mesh.D * (ṽ)
β = weights * (Δxⱼ .^(2 * 1 -1) .* (deriv .^2))
for jjj in 2:n
    deriv .= mesh.D * deriv
    β .+= weights * (Δxⱼ .^(2 * jjj -1) .* (deriv .^2))
end
sumβ = β[i1 .- 1] + β[i1] + β[i1 .+ 1]
γ0 = 0.001
γ1 = 0.998
γ2 = γ0
ϵ = 1e-6
allw = 1 ./ ((ϵ .+ β) .^ 2)
tmpw0 = γ0 ./ allw
tmpw1 = γ1 ./ allw
tmpw2 = γ2 ./ allw
function adjust(x)
    if x==0
        return K
    elseif x==(K+1)
        return 1
    else
        return x
    end
end

lefti = adjust.(i1 .- 1)
righti = adjust.(i1 .+ 1)
sumw = γ0*allw[lefti] + γ1*allw[i1] + γ2*allw[righti]
w0 = γ0*allw[lefti] ./ sumw
w1 = γ1*allw[i1] ./ sumw
w2 = γ2*allw[righti] ./ sumw
w0 = reshape(w0, (1, length(i1)))
w1 = reshape(w1, (1, length(i1)))
w2 = reshape(w2, (1, length(i1)))

plot(x[:,i1], ṽ[:,i1] + v̅[:,i1], label = false, color = :red)
plot!(x[:,i1], w0 .* ṽ[:,lefti] + w1 .* ṽ[:,i1]+ w2 .* ṽ[:,righti]+ v̅[:,i1], label = false, color = :green)

plot(x[:, i2], v[:,i2], color = :blue)
plot!(x[:,i1], w0 .* ṽ[:, lefti] + w1 .* ṽ[:,i1] + w2 .* ṽ[:,righti] + v̅[:,i1], label = false, color = :green)
##
function smoothness_indicator(v, mesh)
    smoothness = zeros(mesh.K)
    Δxⱼ = reshape(mesh.x[end,:]-mesh.x[1,:], (1,mesh.K))
    deriv = mesh.D * (v)
    weights = sum(mesh.M, dims = 1)
    β = weights * (Δxⱼ .^(2 * 1 -1) .* (deriv .^2))
    for jjj in 2:n
        deriv .= mesh.D * deriv
        β .+= weights * (Δxⱼ .^(2 * jjj -1) .* (deriv .^2))
    end
    return β
end

function troubled_cells(v, mesh, cells)
    K = mesh.K
    v̅ = cells * v
    ṽ = v̅[1,:] -  v[1,:]
    ṽ̃ = v[end, :] - v̅[end,:]
    vtb = reshape(v̅[mesh.vmapP] - v̅[mesh.vmapM], (2,80))
    Δ₊v = vtb[end,:]
    Δ₋v = -vtb[1,:]
    troubled_i = collect(1:K)
    minmods = [minmod([ṽ[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
    troubled = minmods .!= ṽ
    minmods = [minmod([ṽ̃[i],Δ₊v[i],Δ₋v[i]]) for i in 1:K]
    troubled2 = minmods .!= ṽ̃
    makeitdouble = troubled .| troubled2
    i1 = troubled_i[makeitdouble]
    return i1, v̅
end

function weno_adjustment(v, mesh, cells)
    β = smoothness_indicator(v, mesh)
    ϵ = 1e-6
    allw = 1 ./ ((ϵ .+ β) .^ 2)
    γ1 = 0.998
    γ0 = (1-γ1)/2
    γ2 = γ0
    i1, v̅ = troubled_cells(v, mesh, cells)
    lefti = adjust.(i1 .- 1)
    righti = adjust.(i1 .+ 1)
    sumw = γ0*allw[lefti] + γ1*allw[i1] + γ2*allw[righti]
    w0 = γ0*allw[lefti] ./ sumw
    w1 = γ1*allw[i1] ./ sumw
    w2 = γ2*allw[righti] ./ sumw
    w0 = reshape(w0, (1, length(i1)))
    w1 = reshape(w1, (1, length(i1)))
    w2 = reshape(w2, (1, length(i1)))
    ṽ  = v-v̅
    adjustedv = w0 .* ṽ[:, lefti] + w1 .* ṽ[:,i1] + w2 .* ṽ[:,righti] + v̅[:,i1]
    v[:, i1] .= adjustedv
    return nothing
end

## WENO Version, wenostart
# Initial Condition
u⁰(x, a, b) = exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2);
u⁰(x,a,b) = sin(2π/(b-a) * x) + 0.5
T = 1.5
inexact = true
cfl = 0.1
K = 80   # Number of elements
n = 8    # Polynomial Order

r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
avg = Diagonal(zeros(n+1))
avg[1] = 1
cells = V * avg * inv(V)

mesh = create_mesh(Ω, elements = K, polynomial_order =  n) # Generate Uniform Periodic Mesh
x = mesh.x
weights = sum(mesh.M, dims = 1)
Δx = x[2] - x[1]
plot(x, u⁰.(x, Ω.a,  Ω.b))
if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end


u0 = @. u⁰(x, Ω.a, Ω.b) # use initial condition for array
α = 1.5/2#1.5/2; # Rusanov parameter
Δu = 0.5 / n# regularization parameter
field_md = DGMetaData(mesh, nothing, nothing); # wrap field metadata
central = DGMetaData(mesh, u0, Rusanov(0.0));  # wrap derivative metadata
rusanov = DGMetaData(mesh, u0, Rusanov(α));    # wrap derivative metadata
u̇ = Data(u0);
u = Field(Data(u0), field_md);
κbase = minimum([1.5e-8, Δx * α ]) # Diffusivity Constant
κ0 = x .* 0 .+ κbase
kmax = Δx * α * 0 /10 + κbase
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
    u̇ == -∂xᴿ(u*u)*0.5  + ∂xᶜ( κ * ∂xᶜ(u)),
]

p = (pde_equation, u, κ, n);

function dg_burgers_weno!(v̇ , v, params, t)
    # unpack parameters
    equations = params[1]
    u = params[2]
    κ = params[3]
    n = params[4]
    equations[1].lhs.data .=  real.(v[1:(n+1),:])
    weno_adjustment(view(v,(n+2):2n+2,:), mesh, cells)
    equations[2].lhs.data .=  real.(v[(n+2):end,:])
    v̇[1:(n+1),:]   .= compute(equations[1].rhs)
    v̇[(n+2):end,:] .= compute(equations[2].rhs)
    return nothing
end

rhs! = dg_burgers_weno!
tspan = (0.0, T)

# Define ODE problem
v = vcat(copy(κ0), copy(u0))
v̇ = (copy(κ0), copy(u0))
ode_problem = (rhs!, v, tspan, p);

using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Euler() # Heun(), RK4, Tsit5, Feagin14(), SSPRK33
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
    v = real.(sol.u[i])[(n+2):end,:]
    i1, v̅ = troubled_cells(v, mesh, cells)
    plt = plot(x, v, label = false, color = :blue, ylims = (-0.6, 1.6))
    plot!(x[:,i1], v[:,i1], label = false, color = :red)
    display(plt)
    sleep(0.05)
end
