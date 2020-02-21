using SpecialFunctions # for gamma function reasons
using LinearAlgebra    # for Guass quadrature

using Revise
using OffsetArrays

"""
unimesh1D(xmin, xmax, K)
# Description
    Generates a uniform 1D mesh
# Arguments
    xmin: smallest value of array
    xmax: largest values of array
    K: number of elements in an array
# Return Values: VX, EtoV
    VX: vertex values | an Array of size K+1
    EtoV: element to node connectivity | a Matrix of size Kx2
# Example
xmin = -1
xmax =  1
K    =  4
VX, EtoV = unimesh1D(xmin, xmax, K)
"""
function unimesh1D(xmin, xmax, K)
    VX = collect(0:K) ./ K .* (xmax - xmin) .+ xmin
    EtoV = Int.(ones(K, 2))
    for i = 1:K
        EtoV[i,1] = Int(i)
        EtoV[i,2] = Int(i+1)
    end
    return VX, EtoV
end

# Mathy aliases
const Γ = gamma

# Coefficients in the Jacobi polynomial recurrence relations.
aᴾ(α, β, n) = 2/(2n+α+β) * √(n * (n+α+β) * (n+α) * (n+β) / (2n+α+β-1) / (2n+α+β+1))
bᴾ(α, β, n) = -(α^2 - β^2) / (2n+α+β) / (2n+α+β+2)

#code checked against the matlab code
"""
jacobi(x, α, β, n)
# Description
- Evaluates the jacobi polynomial at the point x
# Arguments
- `x`: point at which you will evaluate the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order
# Return
-  `y`: the value of the of the Jacobi polynomial
"""
function jacobi(x, α, β, n::Int)
    Pᵅᵝ = n <= 1 ? OffsetArray(zeros(2), 0:1) : OffsetArray(zeros(n+1), 0:n)
    Pᵅᵝ[0] = √(2.0^-(α+β+1) * Γ(α+β+2) / Γ(α+1) / Γ(β+1))
    Pᵅᵝ[1] = Pᵅᵝ[0]/2 * √((α+β+3) / (α+1) / (β+1)) * ((α+β+2)*x + α - β)
    for n′ in 1:n-1
        Pᵅᵝ[n′+1] = ((x - bᴾ(α,β,n′)) * Pᵅᵝ[n′] - aᴾ(α,β,n′) * Pᵅᵝ[n′-1]) / aᴾ(α, β, n′+1)
    end
    return Pᵅᵝ[n]
end

"""
djacobi(x, α, β, n)
# Description
- Evaluates the derivative of the jacobi polynomial at the point x
# Arguments
- `x`: point at which you will evaluate the derivative of the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order
# Return
-  `y`: the derivative of the of the Jacobi polynomial
"""
djacobi(x, α, β, n::Int) = √(n * (n+α+β+1)) * jacobi(x, α+1, β+1, n-1)

"""
vandermonde(x, α, β, N)
# Description
    Return vandermonde matrix of order N at the values x
    Allocates a little bit of memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second parameter for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `v`: vandermonde matrix
# Example
    See LegendreTests.jl
"""
function vandermonde(x, α, β, N)
    # compute first two coefficients
    γ0 = 2^(α + β + 1) * factorial(α) * factorial(β) / ((α + β + 1) * factorial(α + β))
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0

    # create view to assign values
    v = zeros(length(x), N+1)
    v1 = view(v, :, 1)
    @. v1 = 1 / sqrt(γ0)

    # explicitly compute second coefficient
    if N == 0
        return v
    end

    v2 = view(v, :, 2)
    @. v2 = ( (α + β + 2) * x/2 + (α - β)/2) / sqrt(γ1)

    if N == 1
        return v
    end

    aʲ = 2 / (2 + α + β) * sqrt((α+1) * (β+1) / (α + β + 3))

    for i in 3:(N+1)
        # get views for ith, i-1th, and i-2th columns
        vi = view(v, :, i)
        vM1 = view(v, :, i-1)
        vM2 = view(v, :, i-2)

        # compute new a and b values
        h1 = 2 * (i-2) + α + β
        aⁱ = 2 / (h1 + 2) * sqrt((i-1) * (i-1 + α + β) * (i-1 + α) * (i-1 + β) / ((h1 + 1) * (h1 + 3)))
        bⁱ = - (α^2 - β^2) / (h1 * (h1 + 2))

        # compute coefficients for ith column
        @. vi = 1 / aⁱ * (-aʲ * vM2 + (x - bⁱ) * vM1)

        # save a coefficient for next iteration
        aʲ = aⁱ
    end

    return v
end

"""
dvandermonde(x, α, β, N)
# Description
    Return the gradient of the vandermonde matrix of order N at the values x
    Allocates a little bit of memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `vr`: gradient of vandermonde matrix
# Example
    See LegendreTests.jl
"""
function dvandermonde(x, α, β, N)
    # create empty matrix (also handles first set of derivatives)
    vr = zeros(length(x), N+1)

    if N == 0
        return vr
    end

    # set values using vandermonde matrix
    v = vandermonde(x, α+1, β+1, N)
    for i in 1:N
        vi = view(v, :, i)
        vrP1 = view(vr, :, i+1)
        @. vrP1 = sqrt(i * (α + β + i+1)) * vi
    end

    return vr
end

"""
dmatrix(x, α, β, N)
# Description
    Return the differentiation matrix of order N at the values x
    Allocates too much memory
# Arguments
-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include
# Return Values
-   `D`: the differentiation matrix
# Example
    See LegendreTests.jl
"""
function dmatrix(x, α, β, N)
    # calculate vandermonde matrix and grad of vandermonde matrix
    vr = dvandermonde(x, α, β, N)
    v  =  vandermonde(x, α, β, N)

    # calculate values using D = vr * v^-1
    d = vr / v

    return d
end

"""
lift1D(V, y)
for computing fluxes
helps compute a surface integral of a quantity
note that the parentheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return V * (transpose(V) * E)
end

"""
lift1D_v2(V, y)
for computing fluxes
nodal form
helps compute a surface integral of a quantity
note that the parantheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D_v2(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return E
end


"""
jacobiGQ(α, β, N)
# Description
    Guass Quadrature points and weights for the Jacobi Polynomial (α,β)
# Input
α, β: Jacobi polynomial descriptors
N:    order of quadrature points
# Return: x,w
x: quadrature points | array of size N+1
w: quadrature weights | array of size N+1
#Example
α = 0
β = 0
N = 4
x, w = jacobiGQ(α, β, N)
"""
function jacobiGQ(α, β, N)
    N == 0 && return [(α-β) / (α+β+2)], [2]

    # Form symmetric matrix from recurrence.
    dv = OffsetArray(zeros(N+1), 0:N)  # diagonal vector
    ev = OffsetArray(zeros(N+1), 0:N)  # sub/super-diagonal vector

    for n in 0:N
        dv[n] = bᴾ(α, β, n)
        ev[n] = aᴾ(α, β, n)
    end

    # Create full matrix combining the two.
    # Need to pass arrays that are not offset.
    J = SymTridiagonal(dv[0:N], ev[1:N])
    (α + β) ≈ 0 && (J[1, 1] = 0)

    # Compute quadrature points and weights by eigenvalue solve.
    x, V = eigen(J)
    w = @. V[1, :]^2 * 2^(α+β+1) / (α+β+1)
    @. w *= factorial(α) * factorial(β) / factorial(α+β)

    return x, w
end

"""
jacobiGL(α, β, N)
# Description
    Guass Labatto quadrature points for the Jacobi Polynomial (α,β)
    The quadrature weights are computed as well (but not returned)
# Arguments
- `α, β`: Jacobi polynomial descriptors
- `N`:    order of quadrature
# Return: x
- `x`: quadrature points  | array of size N+1
# Examples
```julia-repl
julia> x = jacobiGL(0, 0, 4)
5-element Array{Float64,1}:
 -1.0
 -0.6546536707079759
  4.440892098500626e-16
  0.6546536707079771
  1.0
```
"""
function jacobiGL(α, β, N)
    N == 0 && error("What are you doing? Gauss-Lobatto points only make sense if N >= 1.")
    N == 1 && return [-1, 1]

    x = zeros(N+1)
    x[1], x[N+1] = -1, 1

    x_GQ, _ = jacobiGQ(α+1, β+1, N-2)
    x[2:N] .= x_GQ

    return x
end

# low storage Runge-Kutta coefficients
rk4a = [ 0.0, -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0]
rk4b = [ 1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0]
rk4c = [ 0.0, 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0]

"""
rk_solver!(u̇, u, params, t)
# Description
    time stepping with 4th order runge-kutta
# Arguments
-   `u̇ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters
-   `t`: time to evaluate at
"""
function rk_solver!(rhs!, fields, fluxes, params, dt, Nsteps; auxils = [])
    # Runge-Kutta residual storage
    solutions = []
    for 𝑓 in fields
        ϕᵗ = similar(𝑓.ϕ)
        @. ϕᵗ = 𝑓.ϕ
        push!(solutions, [ϕᵗ])
    end

    # time step loop
    for tstep in 1:Nsteps
        time = dt * tstep
        for iRK in 1:5
            # get numerical solution
            if isempty(auxils)
                rhs!(fields, fluxes, params, time)
            else
                rhs!(fields, fluxes, auxils, params, time)
            end

            # update solutions
            for 𝑓 in fields
                @. 𝑓.r = rk4a[iRK] * 𝑓.r + 𝑓.ϕ̇ * dt
                @. 𝑓.ϕ = rk4b[iRK] * 𝑓.r + 𝑓.ϕ
            end
        end

        for (i,𝑓) in enumerate(fields)
            ϕᵗ = similar(𝑓.ϕ)
            @. ϕᵗ = 𝑓.ϕ
            push!(solutions[i], ϕᵗ)
        end

        if (tstep % 1000) == 0
            println( string(tstep, " / ", Nsteps))
        end
    end

    return solutions
end

"""
rel_error(u,v)
# Description
- calculate the relative error between u and v with respect to v
# Arguments
- `u` : a structure of numbers
- `v` : a structure of numbers
# return
- `relative error`:
"""
function rel_error(u,v)
    return maximum(abs.(u[:] .- v[:])) / maximum(abs.(u[:]))
end

"""
rel_1_error(u,v)
# Description
- calculate the relative error between u and v with respect to v
# Arguments
- `u` : a structure of numbers
- `v` : a structure of numbers
# return
- `relative error`:
"""
function rel_1_error(u,v)
    return sum(abs.(u[:] .- v[:])) / sum(abs.(u[:]))
end



"""
dropϵzeros!(sparseMatrix)
# Description
- Drops machine zeros in sparse matrix
# Arguments
- `!A`: a sparse matrix
# return
- nothing
"""
function dropϵzeros!(A)
    i,j = findnz(A)
    drop_criteria = eps(maximum(abs.(A)))
    for loop in 1:length(i)
        if abs(A[i[loop],j[loop]]) < drop_criteria
            A[i[loop],j[loop]] = 0.0
        end
    end
    dropzeros!(A)
end

"""
dropϵzeros!(sparseMatrix, drop_criteria)
# Description
- Drops machine zeros in sparse matrix
# Arguments
- `A`: a sparse matrix
- `drop_criteria`: criteria for dropping entries
# return
- nothing
"""
function dropϵzeros!(A, drop_criteria)
    i,j = findnz(A)
    for loop in 1:length(i)
        if abs(A[i[loop],j[loop]]) < drop_criteria
            A[i[loop],j[loop]] = 0.0
        end
    end
    dropzeros!(A)
end
