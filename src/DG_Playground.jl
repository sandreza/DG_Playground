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

## Structs
export  Mesh,
        Field1D,
        Gradient,
        Flux,
        Field
## Fluxes / Boundary Conditions
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

## Operators
export  build,
        build_operator,
        ⋅,
        ⊗,
        compute_volume_terms,
        compute_surface_terms
## Basic DG Atoms
export  jacobiGL,
        dmatrix,
        vandermonde,
        lift1D

end # module
