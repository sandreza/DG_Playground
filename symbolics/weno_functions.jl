## WENO FUNCTIONS
# for periodic domains, needs to be fixed for domains with boundaries
function adjust(x, mesh)
    if x==0
        return mesh.K
    elseif x==(mesh.K+1)
        return 1
    else
        return x
    end
end

function minmod(a)
    if prod(sign.(a) .== sign(a[1]))
        return sign.(a[1])*minimum(abs.(a))
    else
        return -0
    end
end
mminmod(a; M = 0.0) = abs.(a[1]) < M ? a[1] : minmod(a)


### Start of Old Version
#=
function smoothness_indicator(v, mesh)
    smoothness = zeros(mesh.K)
    Δxⱼ = reshape(mesh.x[end,:]-mesh.x[1,:], (1,mesh.K))
    deriv = mesh.D * (v)
    weights = sum(mesh.M, dims = 1) ./ 2
    β = weights * (Δxⱼ .^(2 * 1 -1) .* (deriv .^2))
    for jjj in 2:n
        deriv .= mesh.D * deriv
        β .+= weights * (Δxⱼ .^(2 * jjj -1) .* (deriv .^2))  # excessive
    end
    return β
end

function troubled_cells(v, mesh, cells; M = 0.0)
    K = mesh.K
    v̅ = cells * v
    ṽ = v̅[1,:] -  v[1,:]
    ṽ̃ = v[end, :] - v̅[end,:]
    vtb = reshape(v̅[mesh.vmapP] - v̅[mesh.vmapM], (2,K))
    Δ₊v = vtb[end,:]
    Δ₋v = -vtb[1,:]
    troubled_i = collect(1:K)
    minmods = [mminmod([ṽ[i],Δ₊v[i],Δ₋v[i]], M = M) for i in 1:K]
    troubled = minmods .!= ṽ
    minmods = [mminmod([ṽ̃[i],Δ₊v[i],Δ₋v[i]], M = M) for i in 1:K]
    troubled2 = minmods .!= ṽ̃
    makeitdouble = troubled .| troubled2
    i1 = troubled_i[makeitdouble]
    return i1, v̅
end

function weno_adjustment(v, mesh, cells; M = 0)
    β = smoothness_indicator(v, mesh)
    ϵ = 1e-6
    allw = 1 ./ ((ϵ .+ β) .^ 2)
    γ1 = 0.998
    γ0 = (1-γ1)/2
    γ2 = γ0
    i1, v̅ = troubled_cells(v, mesh, cells, M=M)
    lefti = adjust.(i1 .- 1)
    righti = adjust.(i1 .+ 1)
    sumw = γ0*allw[lefti] + γ1*allw[i1] + γ2*allw[righti]
    w0 = γ0*allw[lefti] ./ sumw
    w1 = γ1*allw[i1] ./ sumw
    w2 = γ2*allw[righti] ./ sumw
    w0 = reshape(w0, (1, length(i1)))
    w1 = reshape(w1, (1, length(i1)))
    w2 = reshape(w2, (1, length(i1)))
    ṽ  = v-v̅
    adjustedv = w0 .* ṽ[:, lefti] + w1 .* ṽ[:,i1] + w2 .* ṽ[:,righti] + v̅[:,i1]
    v[:, i1] .= adjustedv
    return nothing
end

function weno_everything(v, mesh, cells; M=0)
    β = smoothness_indicator(v, mesh)
    ϵ = 1e-6
    allw = 1 ./ ((ϵ .+ β) .^ 2)
    γ1 = 0.998
    γ0 = (1-γ1)/2
    γ2 = γ0
    i1 = collect(1:mesh.K)
    v̅ = cells * v
    lefti = adjust.(i1 .- 1)
    righti = adjust.(i1 .+ 1)
    sumw = γ0*allw[lefti] + γ1*allw[i1] + γ2*allw[righti]
    w0 = γ0*allw[lefti] ./ sumw
    w1 = γ1*allw[i1] ./ sumw
    w2 = γ2*allw[righti] ./ sumw
    w0 = reshape(w0, (1, length(i1)))
    w1 = reshape(w1, (1, length(i1)))
    w2 = reshape(w2, (1, length(i1)))
    ṽ  = v-v̅
    adjustedv = w0 .* ṽ[:, lefti] + w1 .* ṽ[:,i1] + w2 .* ṽ[:,righti] + v̅[:,i1]
    v[:, i1] .= adjustedv
    return nothing
end

function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno = true)
    uⁿ = copy(equations[1].lhs.data)
    # weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    equations[1].lhs.data .= u¹
    # weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    equations[1].lhs.data .= u²
    weno_adjustment(u², mesh, cells, M = M) #works?
    #weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end

function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno = true)
    uⁿ = copy(equations[1].lhs.data)
    weno_adjustment(uⁿ, mesh, cells, M = M)
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    equations[1].lhs.data .= u¹
    # weno_adjustment(u¹, mesh, cells, M = M)
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    equations[1].lhs.data .= u²
    # weno_adjustment(u², mesh, cells, M = M) #works?
    #weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end
=#
### End of Old Version

### Start of new version

function smoothness_indicator(v, cellindex, mesh)
    Δxⱼ = mesh.x[end,cellindex]-mesh.x[1,cellindex]
    deriv = mesh.D * v
    weights = sum(mesh.M, dims = 1) ./ 2
    β = weights * (Δxⱼ .^(2 * 1 -1) .* (deriv .^2))
    for jjj in 2:n
        deriv .= mesh.D * deriv
        β .+= weights * (Δxⱼ .^(2 * jjj -1) .* (deriv .^2)) # excessive
    end
    return β
end

function troubled_cells(v, mesh, cells; M = 0.0)
    K = mesh.K
    v̅ = cells * v
    ṽ = v̅[1,:] -  v[1,:]
    ṽ̃ = v[end, :] - v̅[end,:]
    vtb = reshape(v̅[mesh.vmapP] - v̅[mesh.vmapM], (2,K))
    Δ₊v = vtb[end,:]
    Δ₋v = -vtb[1,:]
    troubled_i = collect(1:K)
    minmods = [mminmod([ṽ[i],Δ₊v[i],Δ₋v[i]], M = M) for i in 1:K]
    troubled = minmods .!= ṽ
    minmods = [mminmod([ṽ̃[i],Δ₊v[i],Δ₋v[i]], M = M) for i in 1:K]
    troubled2 = minmods .!= ṽ̃
    makeitdouble = troubled .| troubled2
    i1 = troubled_i[makeitdouble]
    return i1, v̅
end

function extrapolate(v, mesh, cellindex)
    cellindexr = adjust(cellindex+1, mesh)
    cellindexl = adjust(cellindex-1, mesh)
    
    x = mesh.x
    Iⱼ = (x[end, cellindex ] - x[1, cellindex ])
    Iⱼ₊₁ = (x[end, cellindexr ] - x[1, cellindexr ])
    Iⱼ₋₁ = (x[end, cellindexl ] - x[1, cellindexl ])
    r = jacobiGL(0, 0, size(v)[1]-1)
    extrap_r = (r .+ 2) * (Iⱼ/Iⱼ₊₁)
    extrap_l = (r .- 2) * (Iⱼ/Iⱼ₋₁)
    Vr = vandermonde(extrap_r, 0, 0, n)
    Vl = vandermonde(extrap_l, 0, 0, n)
    extrapr = Vr * inv(V)
    extrapl = Vl * inv(V)
    weights = sum(mesh.M, dims = 1) ./ 2
    vr = extrapr * v[:, cellindexl]
    vr = vr # .- weights*vr .+ v̅ 
    vl = extrapl * v[:, cellindexr]
    vl = vl # .- weights*vl .+ v̅ 
    
    # vr = v[:, cellindexl]
    # vl = v[:, cellindexr]
    return vr, vl
end

function reconstruct!(v, vl, vr, cellindex, mesh, weights; γ1 = 0.998, ϵ=1e-6, threshold = 0.0, variational = false)
    vre = (vr[end]+v[end, cellindex])/2
    vle = (vl[1]+v[1, cellindex])/2
    v̅ = weights * v[:, cellindex]
    vl = vl .- weights * vl .+ v̅
    vr = vr .- weights * vr .+ v̅
    # smooth extrapolations before passing in
    # this would matter if w1-w3 were chosen
    # as the "smoothest"
    w1 = 1.0 # w1 = 0.7, w2 = 0.1
    w2 = 0.0
    w3 = 1-w1-w2
    vltmp = (vl .* w1 + vr .* w2) + w3 .* v[:, cellindex] 
    vrtmp = (vl .* w2 + vr .* w1) + w3 .* v[:, cellindex]
    vl = vltmp
    vr = vrtmp
    β₀ =  smoothness_indicator(vl, cellindex, mesh)
    β₁ =  smoothness_indicator(v[:, cellindex], cellindex, mesh)
    β₂ =  smoothness_indicator(vr, cellindex, mesh)

    if variational
        β₀ .+= 0.01*((vre - vl[end])^2 + (vle - vl[1])^2)
        β₁ .+= 0.01*((vre - v[end, cellindex])^2 + (vle - v[1, cellindex])^2) 
        β₂ .+= 0.01*((vre - vr[end])^2 + (vle - vr[1])^2) 
    end
    β = [β₀ β₁ β₂]
    allw = 1 ./ ((ϵ .+ β) .^ 2)
    γ0 = (1-γ1)/2
    γ2 = γ0
    sumw = γ0*allw[1] + γ1*allw[2] + γ2*allw[3]
    w0 = γ0*allw[1] ./ sumw
    w1 = γ1*allw[2] ./ sumw
    w2 = γ2*allw[3] ./ sumw
    if w1 < threshold
        w1 = threshold
        γ0 = γ2 = (1-w1)/2
        w0 = γ0 * allw[1] ./ (allw[1] + allw[3]) 
        w2 = γ2 * allw[3] ./ (allw[1] + allw[3]) 
    end
    v[:, cellindex] .= w0*vl + w1*v[:, cellindex] + w2*vr
    return nothing
end

function weno_adjustment(v, mesh, cells; M = 0, threshold = 0)
    i1, v̅ = troubled_cells(v, mesh, cells; M = M)
    weights = sum(mesh.M, dims = 1) ./ 2
    for i in i1
        vr, vl = extrapolate(v, mesh, i)
        reconstruct!(v, vl, vr, i, mesh, weights, threshold = threshold)
    end
    return nothing
end

#=
function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno = true)
    uⁿ = copy(equations[1].lhs.data)
    weno_adjustment(uⁿ, mesh, cells, M = M)
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    equations[1].lhs.data .= u¹
    weno_adjustment(u¹, mesh, cells, M = M)
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    equations[1].lhs.data .= u²
    weno_adjustment(u², mesh, cells, M = M) #works?
    #weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end
=#

function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true, threshold = 0)
    uⁿ = copy(equations[1].lhs.data)
    equations[1].lhs.data .= uⁿ
    if weno
        weno_adjustment(uⁿ, mesh, cells, M = M,  threshold = threshold)
    end
    equations[1].lhs.data .= uⁿ
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    if weno
        weno_adjustment(u¹, mesh, cells, M = M,  threshold = threshold)
    end
    equations[1].lhs.data .= u¹
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    if weno 
        weno_adjustment(u², mesh, cells, M = M,  threshold = threshold)
    end
    equations[1].lhs.data .= u²
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end

#=
function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true)
    uⁿ = copy(equations[1].lhs.data)
    equations[1].lhs.data .= uⁿ
    if weno     
        weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    end
    
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    equations[1].lhs.data .= u¹
    if weno     
        weno_adjustment(equations[1].lhs.data .= u¹, mesh, cells, M = M)
    end
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    equations[1].lhs.data .= u²
    if weno 
        weno_adjustment(equations[1].lhs.data .= u², mesh, cells, M = M)
    end
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end
=#
#=
function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true, threshold = 0)
    uⁿ = copy(equations[2].lhs.data)
    equations[2].lhs.data .= uⁿ
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end
    u¹ = uⁿ + Δt * compute(equations[2].rhs)

    equations[2].lhs.data .= u¹
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[2].rhs)

    equations[2].lhs.data .= u²
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno 
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end  
    equations[2].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[2].rhs)
    return nothing
end

function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true, threshold = 0)
    uⁿ = copy(equations[2].lhs.data)
    if weno
        weno_adjustment(uⁿ, mesh, cells, M = M,  threshold = threshold)
    end
    equations[2].lhs.data .= uⁿ
    equations[1].lhs.data.data .= compute(equations[1].rhs)

    u¹ = uⁿ + Δt * compute(equations[2].rhs)
    if weno
        weno_adjustment(u¹, mesh, cells, M = M,  threshold = threshold)
    end
    equations[2].lhs.data .= u¹
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[2].rhs)
    if weno
        weno_adjustment(u², mesh, cells, M = M,  threshold = threshold)
    end
    equations[2].lhs.data .= u²
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno 
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end  
    equations[2].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[2].rhs)
    return nothing
end


function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true, threshold = 0)
    uⁿ = copy(equations[2].lhs.data)
    equations[2].lhs.data .= uⁿ
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end
    u¹ = uⁿ + Δt * compute(equations[2].rhs)

    equations[2].lhs.data .= u¹
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[2].rhs)

    equations[2].lhs.data .= u²
    equations[1].lhs.data.data .= compute(equations[1].rhs)
    if weno 
        weno_adjustment(equations[1].lhs.data.data, mesh, cells, M = M,  threshold = threshold)
    end  
    equations[2].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[2].rhs)
    return nothing
end



# incorrect
function ssp3_step!(equations, mesh, cells, Δt; M = 0, weno=true)
    uⁿ = copy(equations[1].lhs.data)
    equations[1].lhs.data .= uⁿ
    if weno     
        weno_adjustment(equations[1].lhs.data, mesh, cells, M = M)
    end
    
    u¹ = uⁿ + Δt * compute(equations[1].rhs)
    equations[1].lhs.data .= u¹
    if weno     
        weno_adjustment(equations[1].lhs.data .= u¹, mesh, cells, M = M)
    end
    u² = 3/4*uⁿ + 1/4*u¹ + (1/4*Δt) * compute(equations[1].rhs)
    equations[1].lhs.data .= u²
    if weno 
        weno_adjustment(equations[1].lhs.data .= u², mesh, cells, M = M)
    end
    equations[1].lhs.data .= 1/3*uⁿ + 2/3*u² + (2/3*Δt) * compute(equations[1].rhs)
    return nothing
end
=#