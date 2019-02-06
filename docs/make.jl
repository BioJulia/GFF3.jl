using Documenter, GFF3

makedocs(
    format = Documenter.HTML(
        edit_branch = "develop"
    ),
    sitename = "GFF3.jl",
    pages = [
        "Home" => "index.md",
    ],
    authors = "Kenta Sato, D. C. Jones, Ben J. Ward, The BioJulia Organisation and other contributors."
)
deploydocs(
    repo = "github.com/BioJulia/GFF3.jl.git",
    deps = nothing,
    make = nothing
)
