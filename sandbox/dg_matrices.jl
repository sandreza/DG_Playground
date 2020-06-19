using DG_Playground

n = 5
α = β = 0.0
# compute Gauss Lobatto grid
r = jacobiGL(α, β, n)

# build differentiation matrix
D = dmatrix(r, α, β, n)

# Vandermonde matrix (the Legendre Transform)
V = vandermonde(r, α, β, n)
# build surface integral terms
lift = lift1D(V)

# build mass matrix and inverse of mass matrix
Mi = V * V'
M = inv(Mi)

S  = M * D
Sᵀ = S'

###
eigvals(M + S * Mi * S')

###
# Checking Legendre modes (prefactors are for normalization)
inv(V) * ( 0 .* r .+ 1)
inv(V) * (r)
inv(V) * (3. * r .* r .- 1 )
inv(V) * (5 .* (r .^ 3) .- 3 .* r)
