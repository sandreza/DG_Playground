
struct test_struct{𝒯}
    a::𝒯
end

a = randn(3)
println(a)
b=test_struct(a)
###

solveAdvection!(v, u, params, 0.0)
@. u += v * dt
solveAdvection!(v, u, params, 0.0)
@. u += v * dt
solveAdvection!(v, u, params, 0.0)
@. u += v * dt
solveAdvection!(v, u, params, 0.0)
@. u += v * dt


###
field_data = sol.u[end]
flux_type = NeglectFlux()
flux_field = Field(field_data, field_bc)
# Flux and state
v = sol.u[end]
Φ = Flux(flux_type, flux_field, v)
p1 = plot(𝒢.x, field_data, legend = false, title = "function")

p2 = plot(𝒢.x, ∇⋅Φ, legend = false, title = "derivative")

# Central
flux_type = Central()
Φ = Flux(flux_type, flux_field, v)
tmp3 = ∇⋅Φ
p3 = plot(𝒢.x, tmp3, legend = false, title = " DG derivative Central")

# Rusanov
flux_type = Rusonov(c)
Φ = Flux(flux_type, flux_field, v)
tmp4 = ∇⋅Φ
p4 = plot(𝒢.x, tmp4, legend = false, title = " DG derivative Rusonov")

plot(p1,p2,p3,p4)
