"""
Core data types for NMR fragment screening analysis
"""
module Types

export Experiment
export ExperimentConfig
export Library
export LibraryFragment
export LibraryCocktail
export Cocktail
export AbstractPeak
export BasicPeak
export LibraryPeak

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
    reference_directory
    experiment_directory
    output_directory
    cocktails
end


struct LibraryFragment
    fragment_id
    smiles
    concentration
    peak_ids
    peak_shifts
    peak_heights
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

abstract type AbstractPeak end

struct LibraryPeak <: AbstractPeak
    id
    refspec
    boundspec
    fragment_id
    cocktail_id
    smiles
    library_shift
end
function Base.show(io::IO, peak::LibraryPeak)
    print(io, "LibraryPeak ", peak.id, ", Cocktail: ", peak.cocktail_id, ", Library Shift (ppm): ", peak.library_shift)
end

struct BasicPeak <: AbstractPeak
    id
    refspec
    boundspec
    fragment_id
    cocktail_id
    smiles
    library_shift
    ref_shift
    bound_shift
    # good::Bool
end
function Base.show(io::IO, peak::BasicPeak)
    print(io, "BasicPeak ", peak.id,
        ", Cocktail: ", peak.cocktail_id,
        ", Library Shift (ppm): ", peak.library_shift,
        ", Ref. Shift (ppm): ", peak.ref_shift,
        ", Bound Shift (ppm): ", peak.bound_shift,
        ", Good: ", peak.good)
end

struct Cocktail
    id
    name
    fragment_ids
    refspec
    boundspec
    peaks::Vector{AbstractPeak}
end

end