module Relaxation

using CairoMakie
using Distances
using FileIO
using GLMakie
using LinearAlgebra
using LsqFit
using Measurements
using MolecularGraph
using NMRTools
using Peaks
using RDKitMinimalLib
using Statistics
using TOML
using UMAP

using ..NMRScreen
using ..Types

export relaxation

# Include other source files
include("types.jl")
include("cocktails.jl")
include("fitting.jl")
include("screeningpeaks.jl")
include("state.jl")
include("save.jl")
include("gui.jl")


function relaxation(config, library, cocktails)
    cocktails = prepcocktails(cocktails)
    state = initialisestate(config, library, cocktails)
    gui!(state)
end

end