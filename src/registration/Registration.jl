module Registration

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

export registration

include("types.jl")
include("nestedobservable.jl")
include("rpm.jl")
include("peaks.jl")
include("cocktails.jl")
include("state.jl")
include("gui.jl")
include("screeningpeaks.jl")


"""
    registration(config, library, cocktails)

Start the registration analysis interface.

# Returns
- `cocktails`: Updated cocktails data structure
- `state`: Updated state data structure (for returning to the GUI)
"""
function registration(config, library, cocktails)
    # 1. prepare registration cocktails
    regcocktails = prepcocktails(cocktails)

    # 2. initialise state for GUI
    state = initialisestate(config, library, regcocktails)

    # 3. show GUI
    registration(state)
end

# start registration from state - for returning to the GUI
function registration(state)
    gui!(state)
    
    # 4. generate updated cocktails from state
    cocktails = recreatecocktails(state)

    return cocktails, state
end

end