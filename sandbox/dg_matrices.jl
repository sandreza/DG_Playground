using DG_Playground

n = 5
α = β = 0.0
# compute Gauss Lobatto grid
r = jacobiGL(α, β, n)

# build differentiation matrix
D = dmatrix(r, α, β, n)

# build surface integral terms
V = vandermonde(r, α, β, n)
lift = lift1D(V)

# build mass matrix and inverse of mass matrix
Mi = V * V'
M = inv(Mi)

S  = M * D
Sᵀ = S'

###
eigvals(M + S * Mi * S')
