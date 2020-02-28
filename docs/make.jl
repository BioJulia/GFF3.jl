using Pkg
using Documenter, GFF3

makedocs(
    format = Documenter.HTML(
        edit_link = "develop"
    ),
    modules = [GFF3],
    sitename = "GFF3.jl",
    pages = [
        "Home" => "index.md",
        "GFF3" => "man/gff3.md",
        "API Reference" => "man/api.md"
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/GFF3.jl.git",
    devbranch = "develop",
    push_preview = true
)
