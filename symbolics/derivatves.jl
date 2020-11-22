include(pwd() * "/symbolics/dynamic_kappa.jl")
##
data = u.data.data
state = data
inexact = false
mesh = create_mesh(Ω, elements = K, polynomial_order =  n)
if inexact
    DM = Diagonal(sum(mesh.M, dims = 1)[:])
    mesh.M .= DM
    mesh.Mi .= inv(DM)
    mesh.lift[:,1] .= mesh.Mi[1,:]
    mesh.lift[:,end] .= mesh.Mi[end,:]
end

method = u.metadata.method
∫dV = compute_volume_terms(data, mesh)
∫dA = compute_surface_terms(mesh, data, state, Rusanov(0))

if inexact
    quad_string = ", inexact quadrature,"
else
    quad_string = ", exact quadrature,"
end
p1 = plot(x, u.data.data, label = false, title = "state" * quad_string * " p="*string(n))
p2 = plot(x, ∫dV, label= false, title = "volume derivative")
p3 = plot(x, ∫dA, label= false, title = "surface derivative")
p4 = plot(x, ∫dV + ∫dA, label = false, title = "total derivative")
plot(p1,p2,p3,p4)

##
plot(x,∫dV)

tt = sum((∫dV - ∫dA) .^2 , dims = 1)
xC = sum(x, dims = 1) / (n+1)
plot(xC[:], tt[:])
##

# projection onto Legendre modes
mesh = create_mesh(Ω, elements = K, polynomial_order =  n)
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
filter = Diagonal(zeros(n+1))
filter[1] = 1
filter[2,2] = 1
filter[3,3] = 0
filter[4,4] = 0
linearfilter = V * filter * inv(V)
plot(x, linearfilter*data,label = false)
plot!(x, data, color = :blue, label = false)

##
# Legendre Mode Amplitudes
lns(x) = log10(abs(x))
p= []
push!(p,scatter(lns.(inv(V) * data[:,1]), label = false, title = "K= " * string(1) ))
# label = "element " * string(1)) 
for i in 2:K
    push!(p,scatter(lns.(inv(V) * data[:,i]), label = false, title = "K= " * string(i) ))
end
ps = plot(p[1:K]...)
pn = plot(x, data,  label = false)
plot(ps,pn)

quant = sum(abs.(spectrum[2:3,:]), dims =1) ./ 5
abs.(spectrum[4,:]) .> quant[:]
abs.(spectrum[5,:]) .> quant[:]
abs.(spectrum[6,:]) .> quant[:]