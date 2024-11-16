"""
Configuration handling for NMR fragment screening
"""


function parseconfig(filename)
    isfile(filename) || throw(SystemError("Configuration file not found: $filename"))
    
    config = TOML.parsefile(filename)
    
    # Validate required sections
    for section in ["setup", "files", "protein", "cocktails"]
        haskey(config, section) || throw(ArgumentError("Missing required section: $section"))
    end
    
    setup = config["setup"]
    files = config["files"]
    protein = config["protein"]
    cocktails = config["cocktails"]
    
    # Required fields in experiment section
    required_setup_fields = ["name", "date", "operator", "spectrometer", "temperature"]
    required_files_fields = ["experiment_directory", "output_directory", "fragment_library"]
    required_protein_fields = ["name", "concentration", "concentration_unit", "buffer"]
    required_cocktail_fields = ["id", "experiment", "reference"]
    
    for field in required_setup_fields
        haskey(setup, field) || throw(ArgumentError("Missing required field in setup section: $field"))
    end
    for field in required_files_fields
        haskey(files, field) || throw(ArgumentError("Missing required field in files section: $field"))
    end
    for field in required_protein_fields
        haskey(protein, field) || throw(ArgumentError("Missing required field in protein section: $field"))
    end
    for cocktail in cocktails
        for field in required_cocktail_fields
            haskey(cocktail, field) || throw(ArgumentError("Missing required field in cocktail section: $field"))
        end
    end

    # get location of config file and use it as working directory
    wd = dirname(filename)
    
    exptconfig = ExperimentConfig(
        setup["name"],
        setup["date"],
        setup["operator"],
        setup["spectrometer"],
        setup["temperature"],
        protein["name"],
        protein["concentration"],
        protein["concentration_unit"],
        protein["buffer"],
        joinpath(wd, files["fragment_library"]),
        joinpath(wd, files["experiment_directory"]),
        joinpath(wd, files["output_directory"]),
        cocktails
    )
        
    # Create output directory if it doesn't exist
    mkpath(exptconfig.output_directory)

    return exptconfig
end
