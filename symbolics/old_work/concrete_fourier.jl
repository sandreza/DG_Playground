include(pwd()*"/symbolics"*"/fourier.jl")

function fourier_nodes(a, b, N)
    return (b-a) .* collect(0:(N-1))/N .+ a
end

function fourier_wavenumbers(a, b, N)
    up = collect(0:1:N-1)
    down = collect(-N:1:-1)
    indices = up
    indices[floor(Int, N/2):end] = down[floor(Int, N/2):end]
    wavenumbers = 2π/(b-a) .* indices
    return wavenumbers
end

function fourier_derivative(y, P, k)
    tmp = copy(y)
    dy = copy(y)
    mul!(tmp, P, y)
    @. tmp *= im * k 
    ldiv!(dy, P, tmp)
    return dy
end