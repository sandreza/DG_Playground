include("advection.jl")
###
field_data = sol.u[end]
state = copy(sol.u[end])
flux_type = NeglectFlux()
state = sol.u[end]

flux_field = Field(field_data, field_bc)
# Flux
Φ = Flux(flux_type, flux_field, state)
p1 = plot(𝒢.x, field_data, legend = false, title = "function")

p2 = plot(𝒢.x, ∇⋅Φ, legend = false, title = "derivative")

# Central
flux_type = Central()
Φ = Flux(flux_type, flux_field, state)
tmp3 = ∇⋅Φ
p3 = plot(𝒢.x, tmp3, legend = false, title = " DG derivative Central")

# Rusanov
flux_type = Rusanov(c)
Φ = Flux(flux_type, flux_field, state)
tmp4 = ∇⋅Φ
p4 = plot(𝒢.x, tmp4, legend = false, title = " DG derivative Rusonov")

plot(p1,p2,p3,p4)
