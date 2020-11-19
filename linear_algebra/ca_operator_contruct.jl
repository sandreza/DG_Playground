using DG_Playground, LinearAlgebra, SparseArrays, Plots
using Printf
"""
ca_operator_constructor(h, Δt, κ¹, κ², L, K, n; μ = 1.0)

# Description
 Constructions a typical DG linear operator for ca (Convective Adjustment). It is of the form A = μ + ∂ᶻ(κ(h) ∂ᶻ), where h is a parameter that denotes the transition from a diffusivity κ² to κ¹, in a domain [0, L]. K corresponds to the number of elements and n corresponds to the polynomial order. By default Neumann boundary conditions are assumed

# Arguments

- `h`: number. Satisfies h ∈ [0, L]. Chooses where the transition of diffusivity occurs
- `γ`: number. Factor that is proportional to time-step size
- `κ¹`: number. diffusivity in x ∈ [h, L], κ¹ > κ²
- `κ²`: number. diffusivity in x ∈ [0, h]
- `L`: number. Domain size
- `K`: number. Number of elements
- `n`: number. polynomial order

# Keyword Argument
- `μ`: number. Used to scale the identity matrix

# Return
- `A`: matrix representing the operator that needs to be inverted the next timestep

# Comment
The size of the matrix is  K(n+1) x K(n+1).
"""
function ca_operator_constructor(h, γ, κ¹, κ², L, K, n; μ = 1.0, inexact = false, freeflux=false, periodic = false, mass_matrix = false)
    function calculate_hyperbolic_flux(x::Number)
        return x
    end
    function calculate_hyperbolic_flux(x::AbstractArray)
        return x
    end
    # Define right hand side of the differential equation
    function convective_adjustment!(u̇, u, params, t)
        # unpack params
        ∇ = params[1]         # Gradient operator
        Φ = params[2]         # flux term
        κ∇Φ = params[3]       # diffusive state
        Φ.state .= u          # update state
        q = ∇⊗Φ               # calculate gradient
        @. κ∇Φ.state = q      # store flux
        tmp =  ∇⋅κ∇Φ          # calculate tendency
        @. u̇ = tmp            # store it
        return nothing
    end
    xmin = 0.0 # left endpoint of domain
    xmax = L   # right endpoint of domain
    𝒢 = Mesh(K, n, xmin, xmax, periodic = periodic) # Generate Mesh
    if inexact
        mesh = 𝒢
        DM = Diagonal(sum(mesh.M, dims = 1)[:])
        mesh.M .= DM
        mesh.Mi .= inv(DM)
        mesh.lift[:,1] .= mesh.Mi[1,:]
        mesh.lift[:,end] .= mesh.Mi[end,:]
    end
    ∇ = Gradient(𝒢)
    a = 0.0
    b = 0.0
    bc = [a b]
    neumann = true
    # Define hyperbolic flux
    α = 0.0 # Rusanov prameter
    flux_type = Rusanov(α)
    if neumann
        field_bc = FreeFlux()
    else
        field_bc = Dirichlet(bc...)
    end
    if freeflux
        field_bc = FreeFlux()
    end
    if periodic
        field_bc = Periodic()
    end
    u = copy(𝒢.x)
    field_data = copy(u)
    flux_field = Field(field_data, field_bc)
    state = copy(u)
    Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

    # Define Diffusive flux
    α = 0.0 # Rusanov parameter
    flux_type = Rusanov(α)
    if neumann
        field_bc = Dirichlet(bc...)
    else
        field_bc = FreeFlux()
    end
    if freeflux
        field_bc = FreeFlux()
    end
    if periodic
        field_bc = Periodic()
    end
    field_data = copy(u)
    flux_field = Field(field_data, field_bc)
    state = copy(u)
    g(x) = x >= h ? κ¹ : κ²
    κ = g.(𝒢.x)
    function calculate_parabolic_flux(x::AbstractArray, κ)
        return κ .* x
    end
    # This function does not matter since bc is zero
    function calculate_parabolic_flux(x::Number)
        return x
    end
    calculate_parabolic_flux(x) = calculate_parabolic_flux(x, κ)
    κ∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

    # Define Diffusion parameters
    params = (∇, Φ, κ∇Φ)
    rhs! = convective_adjustment!

    affine_operator!(x,y) = rhs!(x, y, params, 0.0)
    x = copy(u)
    Ax = copy(u)
    A, b = build_operator(affine_operator!, 𝒢, mass_matrix = mass_matrix)
    L = μ .* I - γ .* A
    return L, κ
end

###
(γ, κ¹, κ², L, K, n) = (1.0, 1.0, 1.0, 1.0, 10, 3)
simple_operator_constructor(h) = ca_operator_constructor(h, γ, κ¹, κ², L, K, n, μ = 0.0)

# Define operators
vector_space_size = K * (n+1)
# for convenience
all_operators = []
# Define operators L0, L1, L2, ..., LK(n+1)
for v in 0:vector_space_size+1
   Llabel = Meta.parse("L" * string(v))
   κlabel = Meta.parse("κ" * string(v))
   @eval $Llabel, $κlabel  = simple_operator_constructor($v / vector_space_size)
   @eval push!(all_operators, $Llabel)
end

###

distances = zeros(vector_space_size + 2 , vector_space_size + 2)
for i in eachindex(all_operators)
    for j in eachindex(all_operators)
        distances[i, j] = norm(all_operators[i]-all_operators[j])
    end
end

###
threshold(L; ϵ = eps(100.0)) = abs(L) > ϵ ? 0 : 1
theme(:juno)
spy(sparse(L2))
heatmap(reverse(threshold.(L2), dims = 1)) 
heatmap(L2)
# eigvals(L2)
diag(L2)

##
xmin = 0.0 # left endpoint of domain
xmax = L   # right endpoint of domain
# 60, (n, K) = (30, 1) , (20, 2), (15, 3), (12, 4), (10, 5)
K = 12
n = 4
G = Mesh(K, n, xmin, xmax) # Generate Mesh
mass_matrix = false
inexact = true
constructor(h) = ca_operator_constructor(h, -1, κ¹, κ², L, K, n, μ = 0.0, inexact = inexact, freeflux = false, periodic = true, mass_matrix = mass_matrix)
Δ = constructor(L/2)[1]
if mass_matrix
    Δ = Symmetric(Δ)
end
heatmap(reverse(threshold.(Δ), dims = 1)) 
λ, vΔ =  eigen(Δ)
Δ² = Δ * Δ

test = cos.(2π/L * G.x)
dxxtest = reshape(Δ * test[:], (n+1, K))
dx4test = reshape(Δ² * test[:], (n+1, K))

plot(G.x, dxxtest, legend = false)
plot(G.x, dx4test, legend = false)
norm(dxxtest + (2π/L)^2 .* test) ./ ((2π/L)^2 )
norm(dx4test - (2π/L)^4 .* test) ./ ((2π/L)^4 )

plots = []
for i in 1:9 
    ind = length(λ)-(i-1)
    eigval = @sprintf("%1.1f", real(λ[ind]))
    p1 = plot(G.x, real.(reshape(vΔ[:,ind], (n+1, K))), legend = false, title = "λ=" * eigval)
    push!(plots, p1)
end
plot(plots...)
##
p1 = plot(G.x, reshape(vΔ[:,ind], (n+1, K)), 
legend = false, title = "Eigenvector with eigenvalue " * eigval, 
ticklabelsize = 0.1 )
