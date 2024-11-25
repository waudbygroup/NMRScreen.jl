function rpm(expected_peaks::Vector{Float64}, 
                observed_peaks::Vector{Float64};
                position_weight::Float64=1.0,
                match_threshold::Float64=0.25,
                T_init::Float64=1.0, 
                T_final::Float64=0.01, 
                alpha::Float64=0.9)
    
    # Convert vectors to Nx2 and Mx2 matrices with uniform intensities
    expected_peaks_matrix = hcat(expected_peaks, ones(length(expected_peaks)))
    observed_peaks_matrix = hcat(observed_peaks, ones(length(observed_peaks)))
    
    # Call the original rpm function with the new matrices
    return rpm(expected_peaks_matrix, observed_peaks_matrix;
               position_weight=position_weight,
               intensity_weight=0.0,  # Intensity differences are irrelevant
               match_threshold=match_threshold,
               T_init=T_init, 
               T_final=T_final, 
               alpha=alpha)
end

"""
    rpm(expected_peaks, observed_peaks; 
        position_weight=1.0,
        intensity_weight=1.0,
        match_threshold=0.5,
        T_init=1.0, 
        T_final=0.01, 
        alpha=0.9)

Match NMR peaks using Gold's Robust Point Matching algorithm,
allowing x-translation and y-scaling.

Parameters:
- `expected_peaks`: Nx2 matrix of (position, intensity) for expected peaks
- `observed_peaks`: Mx2 matrix of (position, intensity) for observed peaks
- `position_weight`: Weight for chemical shift differences in distance calculation
- `intensity_weight`: Weight for intensity differences in distance calculation
- `match_threshold`: Minimum correspondence value to consider peaks matched
- `T_init`: Initial temperature for deterministic annealing
- `T_final`: Final temperature
- `alpha`: Temperature reduction factor

Returns:
- `x_offset`: Optimal offset in x (position) dimension
- `y_scale`: Optimal scaling in y (intensity) dimension
- `matches`: Vector of (expected_idx, observed_idx) pairs
- `unmatched_expected`: Vector of unmatched expected peak indices
- `unmatched_observed`: Vector of unmatched observed peak indices
- `correspondence`: Final correspondence matrix
"""
function rpm(expected_peaks::Matrix{Float64}, 
                observed_peaks::Matrix{Float64};
                position_weight::Float64=1.0,
                intensity_weight::Float64=0.01,
                match_threshold::Float64=0.25,
                T_init::Float64=1.0, 
                T_final::Float64=0.01, 
                alpha::Float64=0.9)
    
    # Get dimensions
    N, M = size(expected_peaks, 1), size(observed_peaks, 1)
    
    # Initialize parameters
    x_offset = 0.0
    y_scale = 1.0
    T = T_init
    
    # Pre-allocate arrays for efficiency
    D = zeros(N, M)
    correspondence = zeros(N, M)
    Y_transformed = similar(observed_peaks)
    
    # Normalize weights to sum to 1
    total_weight = position_weight + intensity_weight
    pos_w = position_weight / total_weight
    int_w = intensity_weight / total_weight
    
    while T > T_final
        # Transform observed points with current parameters
        Y_transformed[:,1] = observed_peaks[:,1] .- x_offset  # Apply x offset
        Y_transformed[:,2] = observed_peaks[:,2] .* y_scale   # Apply y scaling
        
        # Compute weighted distance matrix
        for i in 1:N, j in 1:M
            # Weighted Euclidean distance with separate position and intensity components
            pos_diff = (expected_peaks[i,1] - Y_transformed[j,1])^2
            int_diff = (expected_peaks[i,2] - Y_transformed[j,2])^2
            D[i,j] = sqrt(pos_w * pos_diff + int_w * int_diff)
        end
        
        # Compute correspondence matrix using softassign
        correspondence .= exp.(-D ./ T)
        
        # Sinkhorn normalization
        for _ in 1:4
            # Normalize rows
            row_sums = sum(correspondence, dims=2)
            correspondence ./= (row_sums .+ 1e-10)
            
            # Normalize columns
            col_sums = sum(correspondence, dims=1)
            correspondence ./= (col_sums .+ 1e-10)
        end
        
        # Update transformation parameters
        total_corr = sum(correspondence)
        
        # Update x_offset (weighted by correspondence and position_weight)
        x_offset_num = 0.0
        for i in 1:N, j in 1:M
            if correspondence[i,j] > 1e-3
                x_offset_num += (observed_peaks[j,1] - expected_peaks[i,1]) * correspondence[i,j]
            end
        end
        x_offset = x_offset_num / (total_corr + 1e-10)
        
        # Update y_scale (weighted by correspondence and intensity_weight)
        y_scale_num = 0.0
        for i in 1:N, j in 1:M
            if correspondence[i,j] > 1e-3
                y_scale_num += (observed_peaks[j,2] / expected_peaks[i,2]) * correspondence[i,j]
            end
        end
        y_scale = y_scale_num / (total_corr + 1e-10)
        
        # Reduce temperature
        T *= alpha
    end
    
    # Get final matches using match_threshold
    matches = Tuple{Int,Int}[]
    used_observed = Set{Int}()
    
    # Find best matches based on final correspondence matrix
    for i in 1:N
        best_match = argmax(correspondence[i,:])
        if correspondence[i,best_match] > match_threshold && !(best_match in used_observed)
            push!(matches, (i, best_match))
            push!(used_observed, best_match)
        end
    end
    
    # Find unmatched peaks
    unmatched_expected = [i for i in 1:N if !any(m[1] == i for m in matches)]
    unmatched_observed = [j for j in 1:M if !(j in used_observed)]

    # visualize_matches(expected_peaks, observed_peaks, x_offset, y_scale, matches) |> display
    # println("Press any key to continue...")
    # readline()
    
    return x_offset, y_scale, matches, unmatched_expected, unmatched_observed #, correspondence
end

# """
#     visualize_matches(expected_peaks, observed_peaks, x_offset, y_scale, matches; 
#                      show_unmatched=true)

# Visualize the peak matches using Plots.jl

# Parameters:
# - `show_unmatched`: If true, highlight unmatched peaks in different color
# """
# function visualize_matches(expected_peaks::Matrix{Float64},
#                          observed_peaks::Matrix{Float64},
#                          x_offset::Float64,
#                          y_scale::Float64,
#                          matches::Vector{Tuple{Int,Int}};
#                          show_unmatched::Bool=true)
    
#     # Create plot
#     p = plot(xlabel="Chemical Shift (ppm)", 
#             ylabel="Intensity",
#             legend=true,
#             grid=true)
    
#     # Get matched indices
#     matched_exp = [m[1] for m in matches]
#     matched_obs = [m[2] for m in matches]
    
#     # Plot expected peaks
#     if show_unmatched
#         # Plot unmatched expected peaks in different color
#         unmatched_exp = setdiff(1:size(expected_peaks,1), matched_exp)
#         if !isempty(unmatched_exp)
#             scatter!(expected_peaks[unmatched_exp,1], expected_peaks[unmatched_exp,2],
#                   label="Expected (unmatched)",
#                   color=:purple,
#                   marker=:circle)
#         end
#     end
    
#     # Plot matched expected peaks
#     scatter!(expected_peaks[matched_exp,1], expected_peaks[matched_exp,2],
#           label="Expected",
#           color=:blue,
#           marker=:circle)
    
#     # Transform observed peaks
#     observed_transformed = copy(observed_peaks)
#     # observed_transformed[:,1] .-= x_offset
#     # observed_transformed[:,2] ./= y_scale
    
#     # Plot transformed observed peaks
#     if show_unmatched
#         # Plot unmatched observed peaks in different color
#         unmatched_obs = setdiff(1:size(observed_peaks,1), matched_obs)
#         if !isempty(unmatched_obs)
#             scatter!(observed_transformed[unmatched_obs,1], observed_transformed[unmatched_obs,2],
#                   label="Observed (unmatched)",
#                   color=:orange,
#                   marker=:circle)
#         end
#     end
    
#     scatter!(observed_transformed[matched_obs,1], observed_transformed[matched_obs,2],
#           label="Observed (transformed)",
#           color=:red,
#           marker=:circle)
    
#     # Draw lines between matches
#     for (exp_idx, obs_idx) in matches
#         plot!([expected_peaks[exp_idx,1], observed_transformed[obs_idx,1]],
#               [expected_peaks[exp_idx,2], observed_transformed[obs_idx,2]],
#               color=:green,
#               linestyle=:dash,
#               alpha=0.5,
#               label="")
#     end
    
#     return p
# end

# # Example usage:
# #=
# # Create sample data
# expected = [10.5 1.0; 15.2 0.5; 20.1 0.8]
# observed = [10.7 0.9; 15.4 0.4; 20.3 0.7; 22.8 0.3]

# # Run matching with custom weights and threshold
# x_offset, y_scale, matches, unmatched_exp, unmatched_obs, correspondence = 
#     rpm_nmr(expected, observed,
#             position_weight=2.0,    # Chemical shift differences count twice as much
#             intensity_weight=1.0,
#             match_threshold=0.4)    # More permissive matching

# # Visualize results
# p = visualize_matches(expected, observed, x_offset, y_scale, matches)
# display(p)
# =#