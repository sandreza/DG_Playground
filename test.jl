

struct mystruct{T}
    a::T
end

function mystruct(a::T) where T
    return mystruct{T}(a)
end

a = mystruct(randn(3))

a.a
