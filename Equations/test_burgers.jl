
include("Equation/Domains.jl")

Ω = PeriodicIntervalDomain(a=0, b=2π) 

burgers_model = :(
    σ = ∂ˣ(u);
    u̇ = -∂ˣ(u * u) + ν * ∂ˣ( σ );
)

balance law = σ = ∂ˣ(u) + ∂ˣ(u * u) + ν * ∂ˣ( σ )

burgers_model = :(balance_law)

eval_symbolics(burgers_model)
# bcs = :(
#     (σ⋅n = 0 , ∂Ω(boundary_ids=...) );
#     (u   = 0 , ∂Ω(boundary_ids=...) );
# )

abstract type AbstractProblem end

struct PDEProblem <: AbstractProblem
    equation
    domain
    bcs
    function PDEProblem(; equation = nothing, domain = nothing, bcs = nothing)
        return new(equation, domain, bcs)
    end
end

problem = PDEProblem(equation = burgers_model, domain = Ω)

# Need to design a flexible type for DGSpatialDiscretizations
dg_model = SpatialDiscretization(
    problem, type = DiscontinuousGalerkin(),
)

abstract type AbstractDGType end
struct DGSEMType <: AbstractDGType
    numerical_flux_first_order::FluxType
    numerical_flux_second_order::FluxType
    numerical_flux_gradient::FluxType
end

DGSEM() = DGSEMType(numerical_flux_first_order = CentralFlux(),
                    numerical_flux_second_order = CentralFlux(),
                    numerical_flux_gradient = CentralFlux())

interpret!(dg_model, DGSEM())

spatial = DiscontinuousGalerkin(∂ˣ,
temporal = RK4(),
state = (:u, :σ),
parameters = ( (ν = 1.0) )
) # automatically annotate with central fluxes everywhere
reinterpret!(equation, 1, :(∂ˣ(u * u)), ::Rusanov(0.1)))

grid = Mesh(Ω, elements = h, nodes = LegendreExtrema(N))
x  = grid.x
u⁰ = exp.(-x^2)
ode_problem = InitialValueProblem(u⁰, equation, grid, Δt = 0.1, adaptive = false)
evolve(ode_problem)

##

K = 20     # Number of elements
n = 2      # Polynomial Order
xmin = 0.0 # left endpoint of domain
xmax = 2π  # right endpoint of domain
mesh = Mesh(K, n, xmin, xmax, periodic = true) 
mesh_concrete = mesh   
x = mesh.x
u = @.  exp(-2 * (xmax-xmin) / 3 * (x - (xmax-xmin)/2)^2)

# (4) Annotate Fluxes
# Define hyperbolic flux (associated with diffusion)
α = -0.0 # Rusanov prameter
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = -0.0 # Rusanov parameter
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)

# Define Advective flux
α = -0.5 # Rusanov parameter (negative)
flux_type = Rusanov(α)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
𝒜Φ = Flux(flux_type, flux_field, state, calculate_advective_flux)