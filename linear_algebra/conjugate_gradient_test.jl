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
include("../dg_utils/tmp_poisson.jl")
solution = copy(sol[:])
b = copy(f[:])
γ = 00.0
s∇² -= γ * I
∇²_tmp(x) = s∇² * x
solution = s∇² \ b
G = inv(∇² - γ*I) # Greens function
G_tmp(x) = G * x
###
# Laplacian
norm(∇²_tmp(solution) - b)/norm(b)
x⁰ = randn(length(sol))
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
x⁰ = randn(length(sol))
P_tmp(x) = G * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Perfect Preconditioner ")

###
# Preconditioned Laplacian, Inverse Diagonal
x⁰ = randn(length(sol))
prec = 1.0 ./ diag(s∇²)
P_tmp(x) = prec .* x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Discrete Tridiagonal inverse Laplacian (-1 2 -1)
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = -2
    if i < length(x⁰)
        prec[i+1,i] = 1
        prec[i,i+1] = 1
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
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = ∇²[i,i] - γ
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
x⁰ = randn(length(sol))
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
nn, mm = size(x)
#=
for i in nn:nn:mm*nn-1
    prec[i+1,i] = 0.0
    prec[i,i+1] = prec[i+1,i]
end
=#
#prec[1,1] = - prec[1,2]
#prec[end,end] = - prec[end, end-1]
# this is a tridiagonal matrix thus doing this is a bad idea
prec = inv(prec)
P_tmp(x) = prec * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(sol))
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
# Preconditioned Laplacian, Tridiagonal band
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = ∇²[i,i]
    if i < length(x⁰)
        prec[i+1,i] = sum(∇²[i+1,:]) - prec[i,i]
        prec[i,i+1] = prec[i+1,i]
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
# Preconditioned Laplacian, Tridiagonal entries of Green's function
x⁰ = randn(length(sol))
prec = Tridiagonal(G)

# this is a tridiagonal matrix thus doing this is a bad idea
P_tmp(x) = prec * x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Bad Preconditioner ")

###
# the solution
scatter(x[:], x⁰)
# within an element, n is the number of GLP in an element
scatter(x[1:(n+1)], x⁰[1:(n+1)])
#
###
# target time
chol_s∇² = cholesky(-s∇²)
@btime chol_s∇² \ b
@btime s∇² * b


###
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = -2
    if i < length(x⁰)
        prec[i+1,i] = 1
        prec[i,i+1] = 1
    end
end
# this is a tridiagonal matrix thus doing this is a bad idea
prec = inv(prec)
mtp = inv(prect)
ind = 3
p1 = plot(x[:], G[ind,:])
p2 = plot(x[:], mtp[ind,:] ./ 10, ylims = (minimum(G[ind,:]), maximum(G[ind,:])))
p3 = plot(x[:], prec[ind,:] ./ 10, ylims = (minimum(G[ind,:]), maximum(G[ind,:])))
plot(p1,p2,p3)
