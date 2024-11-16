"""
Core data types for NMR fragment screening analysis
"""

struct Experiment
    config
    library
    cocktails
end

struct ExperimentConfig
    name
    date
    operator
    spectrometer
    temperature
    protein_name
    protein_concentration
    protein_concentration_unit
    buffer
    library_filename
    experiment_directory
    output_directory
    cocktails
end


struct LibraryFragment
    id
    smiles
    concentration
    peaks
end

struct LibraryCocktail
    id
    name
    fragment_ids
end

struct Library
    name
    fragments::Dict{String,LibraryFragment}
    cocktails::Dict{String,LibraryCocktail}
end


struct Cocktail
    id
    name
    fragments
    refspec
    boundspec
    peaks
    missingpeaks
end


"""
    RelaxationResult

Results from relaxation analysis of a peak.

# Fields
- `R2`: R2 relaxation rate (s⁻¹)
- `R2_error`: Error in R2
- `I0`: Initial intensity
- `I0_error`: Error in I0
"""
struct RelaxationResult
    R2
    R2_error
    I0
    I0_error
    t
    y
    ye
    tpred
    ypred
end

"""
    ScreeningResult

Complete analysis results for a single peak.

# Fields
- `fragment_id`: Fragment identifier
- `peak_number`: Peak number within fragment
- `cocktail_id`: Cocktail identifier
- `library_shift`: Expected chemical shift
- `reference_shift`: Measured shift without protein
- `bound_shift`: Measured shift with protein
- `reference_intensity`: Reference peak intensity (-protein)
- `reference_intensity_error`: Error in reference intensity
- `bound_intensity`: bound peak intensity (-protein)
- `bound_intensity_error`: Error in bound intensity
- `reference_relaxation`: Relaxation analysis without protein
- `bound_relaxation`: Relaxation analysis with protein
"""
struct ScreeningPeak
    # NB constructor is in screeningpeaks.jl
    # fragment_id
    refspec
    boundspec
    peak_id
    cocktail_id
    library_shift
    reference_shift
    bound_shift
    reference_intensity
    reference_intensity_error
    bound_intensity
    bound_intensity_error
    reference_relaxation
    bound_relaxation
end
