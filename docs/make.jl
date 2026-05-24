using Pkg

Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
Pkg.instantiate()

using Documenter
using Literate
using RooFitLite

const ROOT = @__DIR__
const LITERATE = joinpath(ROOT, "literate")
const TUTORIALS = joinpath(ROOT, "src", "tutorials")

mkpath(TUTORIALS)

for tutorial in ("roofit_basics.jl", "signal_background.jl")
    Literate.markdown(
        joinpath(LITERATE, tutorial),
        TUTORIALS;
        documenter=true,
        execute=false,
        credit=false,
    )
end

makedocs(;
    sitename="RooFitLite.jl",
    modules=[RooFitLite],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaHEP.github.io/RooFitLite.jl/stable/",
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "RooFit basics" => "tutorials/roofit_basics.md",
            "Signal plus background" => "tutorials/signal_background.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaHEP/RooFitLite.jl.git",
    push_preview=true,
    versions=["dev" => "dev"],
)
