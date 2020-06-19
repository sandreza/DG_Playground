using DG_Playground, LinearAlgebra, SparseArrays, Plots


"""
simple_operator_constructor(h, Δt, κ¹, κ², L, K, n; μ = 1.0)

# Description
 Constructions a typical DG linear operator

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
function simple_operator_constructor(h, γ, κ¹, κ², L, K, n; μ = 1.0)
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
    𝒢 = Mesh(K, n, xmin, xmax) # Generate Mesh
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
    A, b = build_operator(affine_operator!, 𝒢)
    L = μ .* I - γ .* A
    return L, κ
end


###
(γ, κ¹, κ², L, K, n) = (1.0, 10.0, 1.0, 1.0, 10, 3)
simple_operator_constructor(h) = simple_operator_constructor(h, γ, κ¹, κ², L, K, n, μ = 0.0)

# Define operators
vector_space_size = K * (n+1)
# for convenience
all_operators = []
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
        distances[i,j] = norm(all_operators[i]-all_operators[j])
    end
end
###
theme(:juno)
spy(sparse(L2))
eigvals(L11)
diag(L2)
