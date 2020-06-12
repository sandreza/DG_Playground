module DG_Playground

# Using
using   Revise,
        SparseArrays,
        LinearAlgebra,
        SpecialFunctions,
        Revise,
        OffsetArrays

# Includes
include("./HesthavenWarburton/mesh.jl")
include("./Abstractions/field.jl")
include("./Abstractions/data_structures.jl")
include("./Abstractions/compute_surface.jl")

# Exports
export  Mesh,
        Field1D,
        Gradient,
        Flux,
        Field

export  NeglectFlux,
        Central,
        Upwind,
        Rusanov,
        RusanovBC,
        Slider,
        Dirichlet,
        Dirichlet2,
        Inflow,
        Inflow2,
        Outflow,
        Outflow2,
        Periodic,
        NoFlux,
        FreeFlux

export  build,
        ⋅,
        ⊗,
        compute_volume_terms,
        compute_surface_terms


end # module
