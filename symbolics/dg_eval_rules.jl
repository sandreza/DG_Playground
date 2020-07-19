using SparseArrays, BenchmarkTools, Plots
include(pwd() * "/symbolics/abstract_core.jl")
include(pwd() * "/src" * "/HesthavenWarburton" * "/utils.jl")
include(pwd() * "/src" * "/HesthavenWarburton" * "/mesh.jl")

function create_mesh(Ω::IntervalDomain; elements = K, polynomial_order = n)
    return Mesh(elements, polynomial_order, Ω.a, Ω.b, periodic = Ω.periodic)
end

# Define Numerical Flux Type
struct Rusanov{𝒯}
    α::𝒯
end

function compute_volume_terms(data::AbstractArray, mesh::Mesh)
    q = mesh.D * data
    @. q *= mesh.rx
    return q
end

# Periodic Boundary Conditions
function compute_surface_terms(mesh::AbstractMesh, data, state::AbstractArray, method::Rusanov{𝒯}) where {𝒯, 𝒮}
    # first compute numerical fluxes at interface
    diffs = reshape( (data[mesh.vmapM] + data[mesh.vmapP]), (mesh.nFP * mesh.nFaces, mesh.K))
    # Include factor of 2 for the weak-strong form
    @. diffs *= 1.0 / 2.0
    # Extra dissipation for Rusanov
    @. diffs[:] += method.α * mesh.normals[:] .* (state[mesh.vmapM] - state[mesh.vmapP]) / 2.0
    # Now create jump in flux, (Weak-Strong form)
    @. diffs[:] -= data[mesh.vmapM]
    # Compute Lift Operator
    lifted =  mesh.lift * (mesh.fscale .* mesh.normals .* diffs)
    return lifted
end

function dg_derivative(mesh, data, state, method)
    ∫dV = compute_volume_terms(data, mesh)
    ∫dA = compute_surface_terms(mesh, data, state, method)
    return ∫dV .+ ∫dA
end

struct DGMetaData{𝒮, 𝒯, 𝒰} 
    mesh::𝒮
    state::𝒯
    method::𝒰
end

##

# Derivatives
dg_derivative(y::AbstractArray, md) = dg_derivative(md.mesh, y, md.state, md.method)
dg_derivative(y::AbstractData, md) = dg_derivative(md.mesh, y.data, md.state, md.method)
function eval(e::Gradient{𝒯, 𝒰}) where {𝒯, 𝒰 <: DGMetaData}
    return dg_derivative(eval(e.operand), e.metadata)
end