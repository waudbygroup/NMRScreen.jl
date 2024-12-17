function loadcocktails(config, library)
    cocktails = config.cocktails

    map(cocktail -> loadcocktail(config, library, cocktail), cocktails)
end

function loadcocktail(config, library, cocktail)
    # get cocktail info from experimental parameters (i.e. experiment filenames)
    cocktail_id = cocktail["id"]
    refspec_path = joinpath(config.reference_directory, cocktail["reference"])
    boundspec_path = joinpath(config.experiment_directory, cocktail["experiment"])

    # load the NMR data
    refspec = loadnmr(refspec_path)
    boundspec = loadnmr(boundspec_path)
    if NMRTools.metadata(refspec, F1Dim, :window) ≠ NMRTools.metadata(boundspec, F1Dim, :window)
        @warn "Mismatched spectral windows between reference and bound spectra" refspec_path boundspec_path NMRTools.metadata(refspec, F1Dim, :window) NMRTools.metadata(boundspec, F1Dim, :window)
    end

    refspec = setrelaxtimes(refspec, acqus(refspec, :vdlist))
    boundspec = setrelaxtimes(boundspec, acqus(boundspec, :vdlist))

    # get cocktail and fragment info from the library definitions
    librarycocktail = try
        library.cocktails[cocktail_id]
    catch e
        if e isa KeyError
            error("Cocktail $cocktail_id was not found in library. Check the library cocktail definitions.")
        else
            rethrow()  # Re-throw any other unexpected errors
        end
    end
    cocktail_name = librarycocktail.name
    fragment_ids = librarycocktail.fragment_ids

    peaks = Vector{AbstractPeak}()
    for fragment_id in fragment_ids
        fragment = try
            library.fragments[fragment_id]
            catch e
                if e isa KeyError
                    error("Fragment $fragment_id was not found in library. Check the library fragment definitions.")
                else
                    rethrow()  # Re-throw any other unexpected errors
                end
            end
        fragment_smiles = fragment.smiles
        # fragment_concentration = fragment.concentration
        fragment_peak_ids = fragment.peak_ids
        fragment_peak_shifts = fragment.peak_shifts
        # fragment_peak_heights = fragment.peak_heights
        for i=1:length(fragment_peak_ids)
            peak_id = fragment_peak_ids[i]
            library_shift = fragment_peak_shifts[i]
            push!(peaks, LibraryPeak(peak_id,
                refspec,
                boundspec,
                fragment_id,
                cocktail_id,
                fragment_smiles,
                library_shift))
        end
    end

    sort!(peaks, by=x->x.library_shift, rev=true)

    Cocktail(cocktail_id, cocktail_name, fragment_ids, refspec, boundspec, peaks)
end


function Base.show(io::IO, cocktail::Cocktail)
    print(io, "Cocktail ", cocktail.id, ", Fragments: ", length(cocktail.fragment_ids), ", Peaks: ", length(cocktail.peaks))
end
