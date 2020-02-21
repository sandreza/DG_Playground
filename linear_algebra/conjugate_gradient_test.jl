include("conjugate_gradient.jl")


# simple tests (from wikipedia)
A = [4.0 1.0; 1.0 3.0]
b = [1.0; 2.0]
x⁰ = [2.0; 1.0]
solution = A \b

A_tmp(x) = A*x
B = inv(A)
pre_tmp(x) = B*x

x⁰ = [2.0; 1.0]
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 1)
println("the relative error after one iteration is ")
println(norm(x⁰ - solution) / norm(solution))
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 2)
println("the relative error after two iterations  is ")
println(norm(x⁰ - solution) / norm(solution))
x⁰ = [2.0; 1.0]
conjugate_gradient!(A_tmp, x⁰, b, maximum_iterations = 1, P = pre_tmp)
println("the relative error after one iteration with a perfect preconditioner is ")
println(norm(x⁰ - solution) / norm(solution))

###
# More complex text using 1D DG stuff
include("../dg_utils/dg_poisson_operator.jl")
include("../dg_utils/utils.jl")
include("../dg_utils/mesh.jl")
include("../dg_utils/field.jl")
# set polynomial order and number of elements
n = 3
K = 3

# set domain parameters
L    = 2π
xmin = 0.0
xmax = L

# generate mesh variables
𝒢 = Mesh(K, n, xmin, xmax)
x = 𝒢.x  #extract gridpoints

∇², M = constructLaplacian(n = n, K=K, xmin = xmin, xmax = xmax)
println("Almost symmetric so let's make it exact")
println(norm(∇² - ∇²', Inf))

∇² = Symmetric(∇²)
s∇² = sparse(∇²)
solution = @. sin(x[:])
b = s∇² * solution
γ = 00.0
s∇² -= γ * M # for helmholtz
∇²_tmp(x) = s∇² * x
solution = s∇² \ b
G = inv(∇² - γ .* M) # Greens function
G = Symmetric(G)
G_tmp(x) = G * x

###
# Laplacian
norm(∇²_tmp(solution) - b)/norm(b)
x⁰ = randn(length(solution))
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Convergence of Conjugate Gradient with the Laplacian")
###
# Greens function
x⁰ = randn(length(x))
Gb = copy(solution)
r = conjugate_gradient!(G_tmp, x⁰, Gb,  track_residual = true)
println("The relative error is")
println(norm(x⁰-b)/norm(b))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Convergence of Conjugate Gradient with the Green's function")

###
# Preconditioned Laplacian, Perfect
x⁰ = randn(length(solution))
P_tmp(x) = G * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Perfect Preconditioner ")

###
# Preconditioned Laplacian, Inverse Diagonal
x⁰ = randn(length(solution))
prec = 1.0 ./ diag(s∇²)
P_tmp(x) = prec .* x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Discrete Tridiagonal inverse Laplacian (-1 2 -1)
x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = -2
    if i < length(x⁰)
        prec[i+1,i] = 1
        prec[i,i+1] = 1
    end
end
prec = Tridiagonal(prec)
lu_prec = lu(prec)
P_tmp(x) = lu_prec \ x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = ∇²[i,i]
    if i < length(x⁰)
        prec[i+1,i] = ∇²[i+1,i]
        prec[i,i+1] = ∇²[i+1,i]
    end
end
# this is a tridiagonal matrix thus doing this is a bad idea
prec = inv(prec)
P_tmp(x) = prec * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Discrete Tridiagonal inverse Laplacian (-1 2 -1) with Δx
x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
Δx = x[2:end] - x[1:end-1]
Δx = @. Δx * Δx
for i in eachindex(Δx)
    if Δx[i] == 0
        Δx[i] = Δx[i+1]
    end
end
for i in eachindex(x⁰)
    if (i < length(x⁰)) && (i >1)
        prec[i,i] = -1/Δx[i] - 1/Δx[i-1] - γ
    elseif i==1
        prec[i,i] = -1/Δx[i] - 1/Δx[i] - γ
    else
        prec[i,i] = -1/ Δx[i-1] - 1/Δx[i-1] - γ
    end
    if i < length(x⁰)
        prec[i+1,i] = 1 / Δx[i]
        prec[i,i+1] = prec[i+1,i]
    end
end

prec = Tridiagonal(prec)
lu_prec = lu(prec)
P_tmp(x) = lu_prec \ x

r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = ∇²[i,i]
    if i < length(x⁰)
        prec[i+1,i] = ∇²[i+1,i]
        prec[i,i+1] = ∇²[i+1,i]
    end
end
prec = Tridiagonal(prec)
lu_prec = lu(prec)
P_tmp(x) = lu_prec \ x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = ∇²[i,i]
    if i < length(x⁰)
        prec[i+1,i] = sum(∇²[i+1,:]) - prec[i,i]
        prec[i,i+1] = prec[i+1,i]
    end
end

prec = Tridiagonal(prec)
lu_prec = lu(prec)
P_tmp(x) = lu_prec \ x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Tridiagonal entries of Green's function
x⁰ = randn(length(solution))
prec = Tridiagonal(G)

# this is a tridiagonal matrix thus doing this is a bad idea
P_tmp(x) = prec * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")


###
# target time
chol_s∇² = cholesky(-s∇²)
@btime chol_s∇² \ b
@btime s∇² * b

###
s∇² = sparse(∇²)
solution = @. sin(x[:])
b = s∇² * solution
γ = 10.0
s∇² -= γ * M # for helmholtz
∇²_tmp(x) = s∇² * x
solution = s∇² \ b
G = inv(∇² - γ .* M) # Greens function
G_tmp(x) = G * x
# (-M + ∇²)⁻¹ ≈ -M⁻¹ - M⁻¹ ∇² M⁻¹
Mi = inv(γ* M)
P = Mi  + Mi * ∇² * Mi
P = (P + P') ./ 2.0
P_tmp(x) = P * x
x⁰ = randn(length(solution))
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")



###
# Preconditioned Laplacian, Discrete Tridiagonal inverse Laplacian (-1 2 -1) with Δx

function reduce_duplicates(𝒢)
    # n is the polynomial order
    # K is the number of elements
    n = 𝒢.n
    K = 𝒢.K
    nt = (n+1)*K
    ni = nt - (K-1)
    A = zeros(ni,nt)
    di = 𝒢.vmapM[2:2:end-1]
    di2 = 𝒢.vmapP[2:2:end-1]
    ti = collect(1:nt)
    ndi = setdiff(ti,𝒢.vmapM[2:1:end-1])
    rdi  = (n+1):n:(ni-1)
    rndi = setdiff(1:ni, rdi)
    for i in eachindex(rdi)
        A[rdi[i], di[i]] = 0.5
        A[rdi[i], di2[i]] = 0.5
    end
    for i in eachindex(rndi)
        A[rndi[i], ndi[i]] = 1.0
    end
    return A
end

x⁰ = randn(length(solution))
prec = zeros(length(x⁰), length(x⁰))
Δx = x[2:end] - x[1:end-1]
Δx = @. Δx * Δx
for i in eachindex(Δx)
    if Δx[i] == 0
        Δx[i] = Δx[i+1]
    end
end
for i in eachindex(x⁰)
    if (i < length(x⁰)) && (i >1)
        prec[i,i] = -1/Δx[i] - 1/Δx[i-1] - γ
    elseif i==1
        prec[i,i] = -1/Δx[i] - 1/Δx[i] - γ
    else
        prec[i,i] = -1/ Δx[i-1] - 1/Δx[i-1] - γ
    end
    if i < length(x⁰)
        prec[i+1,i] = 1 / Δx[i]
        prec[i,i+1] = prec[i+1,i]
    end
end

prec = Tridiagonal(prec)
lu_prec = lu(prec)
P_tmp(x) = lu_prec \ x

r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")
