using Pkg
Pkg.activate(@__DIR__)  # Activate the docs environment
Pkg.instantiate()       # Install all dependencies

using Documenter
using NMRScreen
ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

makedocs(
    sitename = "NMRScreen.jl",
    modules = [NMRScreen],
    pages=[
        "Home" => "index.md",
    ],
    authors="Chris Waudby",
    format = Documenter.HTML(),
    warnonly = [:missing_docs],
)

deploydocs(;
    repo="github.com/waudbygroup/NMRScreen.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "main"],
)