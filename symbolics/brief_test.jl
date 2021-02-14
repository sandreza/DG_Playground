mesh = create_mesh(Ω, elements = K, polynomial_order =  1) # Generate Uniform Periodic Mesh

mesh.Mi
DM = Diagonal(sum(mesh.M, dims = 1)[:])


mesh.Mi[:,1]

Array(inv(DM)[:,1])

n = 2
r = jacobiGL(0, 0, n)
V = vandermonde(r, 0, 0, n)
inv(V) * (0 .* r .+ 1)
inv(V) * (1 .* r .+ 0)
inv(V) * ( 3. * r .^2 .- 1)

ℱ *  (3. * r .^2 .- 1 + r .+ 3)
