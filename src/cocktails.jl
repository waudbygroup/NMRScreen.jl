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

    libraryfragments = [library.fragments[id] for id in fragment_ids]
    librarypeak_ids, librarypeaks = getlibrarypeaks(libraryfragments)

    refpeaks = detectpeaks(refspec)
    boundpeaks = detectpeaks(boundspec)

    screeningpeaks, missingpeaks = makescreeningpeaks(refspec, boundspec, cocktail_id, librarypeak_ids, librarypeaks, refpeaks, boundpeaks)

    Cocktail(cocktail_id, cocktail_name, fragment_ids, refspec, boundspec, screeningpeaks, missingpeaks)
end


function getlibrarypeaks(libraryfragments)
    # return list of peaks for all fragments
    # generate peak ids based on fragment ids + alphabetical tag (a, b, c...)
    librarypeak_ids = Vector{String}()
    librarypeaks = [] # TODO typing?

    for fragment in libraryfragments
        fragid = fragment.id
        for (i, peak) in enumerate(eachrow(fragment.peaks))
            id = fragid * ('a' + i - 1)
            push!(librarypeak_ids, id)
            push!(librarypeaks, peak)
        end
    end
    librarypeaks = librarypeaks |> stack |> transpose |> collect

    return librarypeak_ids, librarypeaks
end


function detectpeaks(spec, snr_threshold=8)
    x = data(spec, F1Dim)
    y = data(spec[:,1]) / spec[:noise]

    pks = findmaxima(y)     # detect peaks
    peakheights!(pks; min=snr_threshold) # filter list

    peak_shifts = x[pks.indices]
    peak_heights = y[pks.indices]
    
    [peak_shifts peak_heights] # return matrix of peaks
end
