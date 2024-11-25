function loadcocktails(config, library)
    cocktails = config.cocktails

    map(cocktail -> loadcocktail(config, library, cocktail), cocktails)
end

function loadcocktail(config, library, cocktail)
    # get cocktail info from experimental parameters (i.e. experiment filenames)
    cocktail_id = cocktail["id"]
    refspec_path = joinpath(config.experiment_directory, cocktail["reference"])
    boundspec_path = joinpath(config.experiment_directory, cocktail["experiment"])

    # load the NMR data
    refspec = loadnmr(refspec_path)
    boundspec = loadnmr(boundspec_path)
    refspec = setrelaxtimes(refspec, acqus(refspec, :vdlist))
    boundspec = setrelaxtimes(boundspec, acqus(boundspec, :vdlist))

    # get cocktail and fragment info from the library definitions
    librarycocktail = library.cocktails[cocktail_id]
    cocktail_name = librarycocktail.name
    fragment_ids = librarycocktail.fragment_ids

    peaks = Vector{AbstractPeak}()
    for fragment_id in fragment_ids
        fragment = library.fragments[fragment_id]
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
