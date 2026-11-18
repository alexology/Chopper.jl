using Documenter
using Chopper

DocMeta.setdocmeta!(Chopper, :DocTestSetup, :(using Chopper); recursive=true)

makedocs(;
    modules=[Chopper],
    sitename="Chopper.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/Chopper.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Differences and validation" => "validation.md",
        "Chopper.jl and chopper 0.14.1" => "DIFFERENCES.md",
        "Contamination screen" => "CONTAMINANTS.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/Chopper.jl", devbranch="main")
