include(pwd()*"/symbolics" * "/concrete_fourier.jl")
# Test concrete implementation
N = 2^8
a,b = (0, 2π)
x = fourier_nodes(a, b, N)
k = fourier_wavenumbers(a, b, N)
P = plan_fft(x*(1+0im))

∂ˣ(y) = fourier_derivative(y, P, k)

y = @. sin(2π*x)*(1+0im)
z = ∂ˣ(y)
rz = real.(z)
scatter(x, rz, label = "fft derivative" )
plot!(x, 2π.*cos.(2π.*x), label = "exact" )

# Test Abstract implementation
fourier_meta_data = FourierMetaData(N, k, NoFilter, P)
y_fourier = FourierData(y)
field = FourierField(y_fourier, fourier_meta_data)

∂x = Derivative(fourier_meta_data)
fourier_derivative(y::FourierData, P, k) = fourier_derivative(y.data, P, k)

eval((∂x⋅(field * field) ) + field)
eval(∂x⋅(field)) * eval(field)

# now the following will work
eval(∂x⋅(field * field) * field)

second_deriv = eval(∂x⋅(∂x⋅field))
plot(x, real.(second_deriv) ./ (2π)^2)

##
u0 = @. exp(-2 * (b-a) / 3 * (x - (b-a)/2)^2)*(1+0im)

fourier_meta_data = FourierMetaData(N, k, NoFilter, P)
y_fourier = FourierData(u0)
field = FourierField(y_fourier, fourier_meta_data)
plot(x, real.(u0))
u = field
κ = 0.001
fourier_meta_data = FourierMetaData(N, k, NoFilter, P)
y_fourier = FourierData(y .* 0.0 .+ κ)
κ = FourierField(y_fourier, fourier_meta_data)
# Burgers equation rhs
u̇ = -(∂x⋅(u*u)) + κ * ( ∂x⋅(∂x⋅u) )
params = (u̇, u, κ)

function fourier_burgers!(v̇ , v, params, t)
    # unpack params
    u̇ = params[1]           # Gradient operator
    u = params[2]           # flux term
    κ = params[3]           # diffusion constant
    u.data.data .= real.(v)
    v̇ .= eval(u̇)
    return nothing
end

##

rhs! = fourier_burgers!
tspan = (0.0, 20.0)

# Define ODE problem
ode_problem = (rhs!, u0, tspan, params)

##
using DifferentialEquations
prob = ODEProblem(ode_problem...);
# Solve it
ode_method = Heun() # Heun(), RK4, Tsit5
dt = 0.1 / N
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
    plt = plot(x, real.(sol.u[i]), xlims=(a, b), ylims = (-1.1,1.1), marker = 3,    leg = false)
    plot!(x, real.(sol.u[1]), xlims = (a, b), ylims = (-1.1,1.1), color = "red", leg = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(plt)
    sleep(0.1)
end