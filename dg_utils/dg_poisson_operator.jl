# poisson solver functions in a Field1D framework
include("field.jl")

"""
solvePoisson!(u̇, u, params, t)


# Description

    Evaluate the right hand side for poisson's equation

# Example

K = 2^2 # number of elements
n = 2^2-1 # polynomial order
println("The degrees of freedom are ")
println( (n+1)*K)

# domain parameters
xmin = 0.0
xmax = 2π

par_i = Field1D(K, n, xmin, xmax)
par_e = external_params(1.0, 1.0)
periodic = false
params = (par_i, par_e, periodic)

x = par_i.x
u = par_i.u

@. u = sin(par_i.x) # initial condition
u̇ = par_i.u̇

@btime solvePoisson!(u̇, u, params, t)
scatter!(x,u, leg = false)

"""
function solvePoisson!(u̇, u, params, t)
    # unpack params
    𝒢 = params[1] # internal parameters
    ι = params[2]
    ε = params[3] # external parameters
    periodic = params[4] #case parameter
    q = params[5]  #temporary arrray for allocation, same size as u
    dq = params[6] #temporary array for allocation, same size as dq
    τ = params[7]   #penalty parameter

    # Form field differences at faces
    diffs = reshape( (u[𝒢.vmapM] - u[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    #@. ι.flux = 1//2 * diffs * (ε.v * 𝒢.normals - (1 - ε.α) * abs(ε.v * 𝒢.normals))
    @. ι.flux =  diffs / 2

    # Inflow and Outflow boundary conditions
    if !periodic
        uin  = -u[𝒢.vmapI]
        uout = -u[𝒢.vmapO]
        ι.flux[𝒢.mapI]  =  @. (u[𝒢.vmapI] - uin) / 2
        ι.flux[𝒢.mapO]  =  @. (u[𝒢.vmapO] - uout) / 2
    end

    # rhs of the semi-discrete PDE, ∂ᵗu = ∂ˣq, ∂ˣu  = q
    #first solve for q
    mul!(q, 𝒢.D, u)
    @. q *= 𝒢.rx
    lift = 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* ι.flux )
    @. q -= lift
    # Form field differences at faces for q
    diffs = reshape( (q[𝒢.vmapM] - q[𝒢.vmapP]), (𝒢.nFP * 𝒢.nFaces, 𝒢.K ))
    #@. dq = 1//2 * diffs * (ε.v * 𝒢.normals - (1 - ε.α) * abs(ε.v * 𝒢.normals))
    @. dq = 0 #reset dq
    @. dq = diffs
    #impose neumann boundary conditions for q
    #=
    if !periodic
        qin  = q[𝒢.vmapI]
        qout = q[𝒢.vmapO]
        dq[𝒢.mapI]  =  @. (q[𝒢.vmapI] - qin) / 2
        dq[𝒢.mapO]  =  @. (q[𝒢.vmapO] - qout) / 2
    end
    =#
    #modify with τ
    fluxq = @. (dq / 2 + τ * 𝒢.normals * ι.flux)
    # solve for u̇
    mul!(u̇, 𝒢.D, q)
    @. u̇ *=  𝒢.rx
    lift = 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* fluxq )
    @. u̇ -= lift
    tmp =  𝒢.M * u̇ #multiply by mass matrix
    @. u̇ = tmp / 𝒢.rx
    return nothing

end


#builds the matrix (one column at a time)
function constructLaplacian(𝒢, periodic, τ)
    L = zeros(length(𝒢.x), length(𝒢.x))
    ι = Field1D(𝒢)
    # set external parameters
    ϰ = 1.0   # diffusivity constant, doesnt actually enter in for now
    α = 1.0 # 1 is central flux, 0 is upwind, doesnt actually enter in for now
    ε = external_params(ϰ, α)

    @. ι.u = 0.0
    q = similar(ι.u)
    dq = similar(ι.flux)

    params = (𝒢, ι, ε, periodic, q, dq, τ)
    # construct it column by column
    for i in 1:length(𝒢.x)
        ι.u[i] = 1.0
        solvePoisson!(ι.u̇, ι.u, params, 0)
        @. L[:,i] = ι.u̇[:]
        ι.u[i] = 0.0
    end
    return L
end


"""
constructLaplacian(; xmin = 0.0, xmax = 2π, n = 3, k = 3, periodic = false, τ = 1.0, verbose = false)

# Description
- Constructs the DG Laplacian with uniformly spaced elements and arbitrary polynomial order. The boundary conditions are assumed to be homogeneous and dirichlet. Eventually flags will be put in to impose different boundary conditions, but the only one that can be done now is periodic or homogeneous and dirichlet. It also returns the

# Keyword arguments
- 'xmin': (number), Leftmost point in the domain
- 'xmax': (number), Rightmost point in the domain
- 'n': (integer), Polynomial order
- 'K': (integer), number of elements
- 'periodic': boolean, imposes periodic boundary conditions
- 'τ': (number), penalty parameter
- 'verbose': (boolean), just states a few things

# Return

- L: (array), matrix form of the Laplacian in weak form
- M: (array), Technically M^{-1}L is the second derivative, this term is necessary to remove from the second derivative so that the operator is symmetric.

"""
function constructLaplacian(; xmin = 0.0, xmax = 2π, n = 3, K = 3, periodic = false, τ = 1.0, verbose = false)
    if verbose
        println("We have " * string((n+1)*K) * " degrees of freedom")
        println("The domain is from " * string(xmin) * " to " * string(xmax))
        println("The value of the penalty parameter is " * string(τ))
    end
    # generate mesh variables
    𝒢 = Mesh(K, n, xmin, xmax)

    L = zeros(length(𝒢.x), length(𝒢.x))
    ι = Field1D(𝒢)
    # set external parameters
    ϰ = 1.0   # diffusivity constant, doesnt actually enter in for now
    α = 1.0 # 1 is central flux, 0 is upwind, doesnt actually enter in for now
    ε = external_params(ϰ, α)

    @. ι.u = 0.0
    q = similar(ι.u)
    dq = similar(ι.flux)

    params = (𝒢, ι, ε, periodic, q, dq, τ)
    # construct it column by column
    for i in 1:length(𝒢.x)
        ι.u[i] = 1.0
        solvePoisson!(ι.u̇, ι.u, params, 0)
        @. L[:,i] = ι.u̇[:]
        ι.u[i] = 0.0
    end
    tmp = zeros(typeof(1.0),(n+1)*K, (n+1)*K)
    for i in 1:K
        tmpi = (i-1)*(n+1)+1:i*(n+1)
        tmp[tmpi,tmpi] .= 𝒢.M
    end
    d = Diagonal((1.0 ./ 𝒢.rx[:])')
    tmp = tmp
    return L, d*tmp
end
