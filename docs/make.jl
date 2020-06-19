using Documenter
using DG_Playground

makedocs(
    sitename = "DG_Playground",
    format = Documenter.HTML(),
    pages = [
    "Home" => "index.md",
    "Discontinuous Galerkin" => "dg_methods.md",
    "Convective Adjustment" => "convective_adjustment.md",
    "Function Index" => "function_index.md",
    ],
    modules = [DG_Playground]
)

deploydocs(repo = "github.com/sandreza/DG_Playground.git")
