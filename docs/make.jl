using Documenter
using DG_Playground


dg_methods = Any[
    "Home" => "dg_methods.md",
    "Single Element" => "dg_single_element.md",
    "Boundary Conditions" => "boundary_conditions.md",
    "Variational Crimes" => "inexact_quadrature.md"
]

physics = Any[
    "Home" => "physics.md",
    "Convective Adjustment" => "convective_adjustment.md",
]

makedocs(
    sitename = "DG_Playground",
    format = Documenter.HTML(collapselevel = 1),
    pages = [
    "Home" => "index.md",
    "Discontinuous Galerkin" => dg_methods,
    "Physics" => physics,
    "Function Index" => "function_index.md",
    ],
    modules = [DG_Playground]
)

deploydocs(repo = "github.com/sandreza/DG_Playground.git")
