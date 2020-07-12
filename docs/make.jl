using Documenter, WignerDMatrices

DocMeta.setdocmeta!(WignerDMatrices, :DocTestSetup, :(using WignerDMatrices); recursive=true)

makedocs(;sitename="WignerDMatrices",
	modules=[WignerDMatrices],
    authors="Jishnu Bhattacharya",
    repo="https://github.com/jishnub/WignerDMatrices.jl/blob/{commit}{path}#L{line}",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages=[
        "Reference" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jishnub/WignerDMatrices.jl",
)