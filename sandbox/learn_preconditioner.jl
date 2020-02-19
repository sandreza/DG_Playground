# using Optim, Revise
include("../mcmc_utils/mcmc.jl")

using Random
m = length(x)
prec = zeros(m, m)
x⁰ = randn(length(sol))
prec = zeros(length(x⁰), length(x⁰))
for i in eachindex(x⁰)
    prec[i,i] = -2
    if i < length(x⁰)
        prec[i+1,i] = 1
        prec[i,i+1] = 1
    end
end
# mat = Tridiagonal(prec)

function closure_loss_for_greens(y::AbstractArray, ∇²::AbstractArray, m::Number; max_it = m)
    prect = zeros((m,m))
    for i in 1:m
            prect[i,i] = y[i]
        if i < m
            prect[i+1,i] = y[i+m]
            prect[i,i+1] = prect[i+1,i]
        end
    end
    #prect = Tridiagonal(prect)
    #λ⁻ = eigvals(prect \ ∇²) .- 1
    #return mean(abs.(λ⁻))
    #P_tmp(x) = prect \ x
    #Random.seed!(1234)
    #x⁰ = randn(m)
    #r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp, tolerance = eps(10^5.0))
    #return length(r)
    # norm(inv(prect) - inv(∇²))
    prect = Tridiagonal(prect)
    # return norm(inv(prect) - inv(∇²))
    lu_prect = lu(prect)
    P_tmp(x) = lu_prect \ x
    Random.seed!(1234)
    x⁰ = randn(m)
    r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp, maximum_iterations = max_it)
    return r[end]
end

function closure_loss_for_greens_2(y::AbstractArray, ∇²::AbstractArray, m::Number; row_info = zeros(m), max_it = m)
    prect = zeros((m,m))
    for i in 1:(m-1)
        prect[i+1,i] = y[i]
        prect[i,i+1] = prect[i+1,i]
    end
    for i in 2:(m-1)
        prect[i,i] = - prect[i,i+1] - prect[i,i-1] + row_info[i]
    end
    prect[1,1] = -2 * prect[2,1] + row_info[1]
    prect[m,m] = -2 * prect[m,m-1] + row_info[m]
    prect = Tridiagonal(prect)
    # return norm(inv(prect) - inv(∇²))
    lu_prect = lu(prect)
    P_tmp(x) = lu_prect \ x
    Random.seed!(1234)
    x⁰ = randn(m)
    r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp, maximum_iterations = max_it)
    return r[end]
end

ttt = sum(∇², dims = 2)[:] .* 0.0
loss(x) = closure_loss_for_greens(x, ∇², m, max_it = 10)
loss_2(x) = closure_loss_for_greens_2(x, ∇², m, row_info = ttt, max_it = 10)

starting_value = randn(m + m - 1)
for i in 1:m
    starting_value[i] = prec[i,i]
end
for i in 1:m-1
    starting_value[i+m] = prec[i+1,i]
end
y = starting_value
y2 = starting_value[m+1:end]
loss(y)
loss_2(y2)

tmpc, sig = optimize_and_estimate_proposal(y, loss, nt = 100000, restart = 3, proposal = [], scale = 0.01, rescale = true, freq = 1000, verbose = true)
loss_2(y2)
loss_2(tmpc)


prect = zeros((m,m))

for i in 1:m
        prect[i,i] = tmpc[i]
    if i < m
        prect[i+1,i] = tmpc[i+m]
        prect[i,i+1] = prect[i+1,i]
    end
end

#=
for i in 1:(m-1)
    prect[i+1,i] = tmpc[i]
    prect[i,i+1] = prect[i+1,i]
end
for i in 2:(m-1)
    prect[i,i] = - prect[i,i+1] - prect[i,i-1] + ttt[i]
end
prect[1,1] = -2 * prect[2,1] + ttt[1]
prect[m,m] = -2 * prect[m,m-1] + ttt[m]
=#
prect = Tridiagonal(prect)
#display(prec)
#display(prect)
loss(y)
loss(tmpc)
###
x⁰ = randn(m)
P_tmp(x) = prect \ x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Optimized Preconditioner ")

x⁰ = randn(m)
P_tmp(x) = x
r = conjugate_gradient!(∇²_tmp, x⁰, b, track_residual = true, P = P_tmp)
println("The relative error is")
println(norm(x⁰-solution)/norm(solution))
scatter!(log.(r)/log(10), ylabel = "log10 residual norm", xlabel = "iterations", title = "Laplacian with Optimized Preconditioner ")
