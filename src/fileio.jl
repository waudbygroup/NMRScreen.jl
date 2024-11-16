"""
Functions for reading and writing NMR screening data
"""


"""
    write_results(filename::String, results::Vector{ScreeningResult})

Write screening results to a CSV file.

# Arguments
- `filename::String`: Output file path
- `results::Vector{ScreeningResult}`: Vector of analysis results

# Returns
Nothing, but creates a CSV file with columns:
- fragment_id: Fragment identifier
- peak_number: Peak number within fragment
- cocktail_id: Cocktail identifier
- library_shift: Expected chemical shift
- reference_shift: Measured shift without protein
- bound_shift: Measured shift with protein
- chemical_shift_perturbation: Change in chemical shift
- intensity_ratio: Intensity ratio (+protein/-protein)
- intensity_ratio_error: Error in intensity ratio
- reference_R2: R2 without protein
- reference_R2_error: Error in reference R2
- bound_R2: R2 with protein
- bound_R2_error: Error in bound R2
- DeltaR2: Change in R2
- DeltaR2_error: Error in ΔR2
"""
function write_results(filename::String, results::Vector{ScreeningResult})
    df = DataFrame(
        fragment_id = String[],
        peak_number = Int[],
        cocktail_id = String[],
        library_shift = Float64[],
        reference_shift = Float64[],
        bound_shift = Float64[],
        chemical_shift_perturbation = Float64[],
        intensity_ratio = Float64[],
        intensity_ratio_error = Float64[],
        reference_R2 = Float64[],
        reference_R2_error = Float64[],
        bound_R2 = Float64[],
        bound_R2_error = Float64[],
        DeltaR2 = Float64[],
        DeltaR2_error = Float64[]
    )
    
    for result in results
        # Calculate derived parameters
        csp = result.bound_shift - result.reference_shift
        delta_R2 = result.bound_relaxation.R2 - result.reference_relaxation.R2
        delta_R2_error = sqrt(result.bound_relaxation.R2_error^2 + 
                            result.reference_relaxation.R2_error^2)
        
        push!(df, (
            result.fragment_id,
            result.peak_number,
            result.cocktail_id,
            result.library_shift,
            result.reference_shift,
            result.bound_shift,
            csp,
            result.intensity_ratio,
            result.intensity_ratio_error,
            result.reference_relaxation.R2,
            result.reference_relaxation.R2_error,
            result.bound_relaxation.R2,
            result.bound_relaxation.R2_error,
            delta_R2,
            delta_R2_error
        ))
    end
    
    CSV.write(filename, df)
    return nothing
end
