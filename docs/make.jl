using Documenter
using DG_Playground

makedocs(
    sitename = "DG_Playground",
    format = Documenter.HTML(),
    pages = [
    "Home" => "index.md",
    ],
    modules = [DG_Playground]
)

deploydocs(repo = "github.com/sandreza/DG_Playground.git")
