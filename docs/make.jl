using Documenter, TrigPolys

makedocs(
    sitename="TrigPolys",

    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),

    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,

    pages = [
        "Introduction" => "index.md",
    ]
)

deploydocs(
    repo   = "github.com/JuliaAlgebra/TrigPolys.jl.git",

    push_preview = true,
)
