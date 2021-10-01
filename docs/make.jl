using Documenter, TrigPolys

makedocs(
    sitename="TrigPolys",

    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),

    pages = [
        "Introduction" => "index.md",
    ]
)

deploydocs(
    repo   = "github.com/yuanchenyang/TrigPolys.jl.git",
)
