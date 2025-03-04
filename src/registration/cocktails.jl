function prepcocktails(cocktails)
    # detect peaks in spectra and assign fragment peaks
    map(cocktails) do cocktail
        # 1.1. detect peaks
        ids = [p.id for p in cocktail.peaks]
        fragment_ids = [p.fragment_id for p in cocktail.peaks]
        smiles = [p.smiles for p in cocktail.peaks]
        xlibrary = [p.library_shift for p in cocktail.peaks]
        xref0 = detectpeaks(cocktail.refspec)
        xbound0 = detectpeaks(cocktail.boundspec)
        
        # 1.2. assign fragment peaks
        xref, xbound, good = process_peaks(ids, xlibrary, xref0, xbound0)

        # apply offset to bound spectrum
        xoff1, xoff2 = getoffsets(xlibrary, xref, xbound, good)
        # @info "Offsets: $xoff1, $xoff2"
        xref .-= xoff1
        xbound .-= xoff2
        newrefspec = add_offset(cocktail.refspec, F1Dim, -xoff1)
        newboundspec = add_offset(cocktail.boundspec, F1Dim, -xoff2)

        RegistrationCocktail(cocktail.id, cocktail.name, fragment_ids, smiles,
            newrefspec, newboundspec,
            # cocktail.refspec, cocktail.boundspec,
            ids, xlibrary,
            Observable(xref), Observable(xbound), Observable(good))
    end
end


function process_peaks(ids, xlibrary, xref, xbound)
    n = length(ids)
    @assert length(xlibrary) == n "xlibrary must match length of ids"
    
    # Initialize output arrays
    xref_final = copy(xlibrary)         # Default to library positions
    xbound_final = copy(xlibrary)       # Default to library positions
    good = fill(false, length(ids))     # Default all peaks to not good
    
    # Get matches between reference and library peaks
    xoff1, _, ref_lib_matches = rpm(xlibrary, xref)
    # Get matches between bound and reference peaks
    xoff2, _, bound_ref_matches = rpm(xref, xbound)

    xref_final = xref_final .+ xoff1
    xbound_final = xbound_final .+ xoff1 .+ xoff2
    
    # Create lookup dictionary for bound->ref matches
    ref_to_bound = Dict(exp_idx => obs_idx for (exp_idx, obs_idx) in bound_ref_matches)
    
    # Process each reference-library match
    for (lib_idx, ref_idx) in ref_lib_matches
        xref_final[lib_idx] = xref[ref_idx]  # Assign reference peak position
        
        # Check if this reference peak has a matching bound peak
        if haskey(ref_to_bound, ref_idx)
            bound_idx = ref_to_bound[ref_idx]
            xbound_final[lib_idx] = xbound[bound_idx]  # Assign bound peak position
            good[lib_idx] = true  # Mark as good since we have all three matches
        else
            xbound_final[lib_idx] = xref[ref_idx]  # Use reference position for bound
            good[lib_idx] = false  # Mark as not good due to missing bound peak
        end
    end
    
    # Peaks without reference matches already have library positions from initialization
    # and are marked as not good by default
    
    # return the final peak positions and good/bad flags, along with the offset for the bound->ref matches
    return xref_final, xbound_final, good
end