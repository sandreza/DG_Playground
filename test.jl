
struct test_struct{𝒯}
    a::𝒯
end

a = randn(3)
println(a)
b=test_struct(a)
###
v = copy(u)
advection!(v, u, params, 0.0)
@. u += v * dt
advection!(v, u, params, 0.0)
@. u += v * dt
advection!(v, u, params, 0.0)
@. u += v * dt
advection!(v, u, params, 0.0)
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


####
α = 1.0 # Rusanov prameter
flux_type = Rusanov(c)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
Φ = Flux(flux_type, flux_field, state, calculate_hyperbolic_flux)

# Define Diffusive flux
α = 1.0 # Rusanov prameter
flux_type = Rusanov(c)
field_bc = Periodic()
field_data = copy(u)
flux_field = Field(field_data, field_bc)
state = copy(u)
∇Φ = Flux(flux_type, flux_field, state, calculate_parabolic_flux)
params = (∇, Φ, ∇ᴰ, ∇Φ)
# unpack params
∇ = params[1]           # Gradient operator
Φ = params[2]           # flux term
∇ᴰ = params[3]          # Gradient operator
∇Φ = params[4]          # Diffusive state
∇Φ.state .= u           # update state
q = ∇ᴰ⋅Φ                # calculate gradient
Φ.state .= q            # store gradient
tmp =  ∇⋅Φ              # calculate (negative) tendency
###
plot(𝒢.x, u)
plot(𝒢.x, q)
plot(𝒢.x, tmp)
