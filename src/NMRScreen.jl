module NMRScreen

using GLMakie
# using Plots
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using Peaks
using TOML

export screen

# Include other source files
include("types.jl")
include("config.jl")
include("library.jl")
include("cocktails.jl")
include("registration.jl")
include("relaxation.jl")
include("screeningpeaks.jl")
include("state.jl")
# include("fileio.jl")
# include("analysis.jl")
include("gui.jl")

"""
    screen()

Entry point for the NMR fragment screening analysis application.

# Arguments
- `experiment_toml`: Path to the experiment configuration TOML file

# Example
```julia
screen("experiment.toml")
```
"""
function screen(experiment_toml::String)
    config = parseconfig(experiment_toml)   # Load configuration
    library = parselibrary(config.library_filename)    # Load library definitions (fragments + cocktails)
    cocktails = loadcocktails(config, library)         # Load cocktail experimental data

    state = initialisestate(config, library, cocktails)
    # Initialize GUI
    fig = gui(state)
    
    # Display and return figure for interactive use
    display(fig)
    return fig, state
    # wait(display(figure)) # wait to execute the following code until the window is closed
    # GLMakie.closeall()
end

"""
Entry point for command line usage.
Following Julia 1.11 convention for standalone applications.
"""
function (@main)(args)
    try
        # Check for command line arguments
        if length(args) != 1
            println("Usage: nmrscreen <experiment_toml>")
            return 1
        end
        
        # Run main function with provided configuration
        screen(args[1])
        
        # Keep the application running
        GLMakie.wait_for_close()
        
        return 0
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
end

end # module

