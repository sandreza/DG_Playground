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

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
