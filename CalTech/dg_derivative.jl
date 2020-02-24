include("advection.jl")
###
field_data = sol.u[end]
state = copy(sol.u[end])
flux_type = NeglectFlux()
state = sol.u[end]

flux_function(x) = x

flux_field = Field(field_data, field_bc)
# Flux
Φ = Flux(flux_type, flux_field, state, flux_function)
p1 = plot(𝒢.x, field_data, legend = false, title = "function", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)

p2 = plot(𝒢.x, ∇⋅Φ, legend = false, title = "derivative", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)

# Central
flux_type = Central()
Φ = Flux(flux_type, flux_field, state, flux_function)
tmp3 = ∇⋅Φ
p3 = plot(𝒢.x, tmp3, legend = false, title = " DG derivative Central", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)

# Rusanov
flux_type = Rusanov(c)
Φ = Flux(flux_type, flux_field, state, flux_function)
tmp4 = ∇⋅Φ
p4 = plot(𝒢.x, tmp4, legend = false, title = " DG derivative Rusonov", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)

plot(p1,p2,p3,p4)
