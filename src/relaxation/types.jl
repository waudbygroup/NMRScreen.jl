"""
    RelaxationResult

Results from relaxation analysis of a peak.
"""
struct RelaxationResult
    R2
    R2_error
    I0
    I0_error
    ty      # vector of points (t,y)
    tye     # vector of tuples (t,y,yerr)
    typred  # vector of points (t,y) for fitted curve
    reducedchi2
end

"""
    ScreeningResult

Complete analysis results for a single peak.
"""
struct ScreeningPeak <: AbstractPeak
    # NB constructor is in screeningpeaks.jl
    # fragment_id
    refspec
    boundspec
    peak_id
    fragment_id
    smiles
    cocktail_id
    umap_x
    umap_y
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
