function parselibrary(config)
    library_filename = config.library_filename
    config = TOML.parsefile(library_filename)
    wd = dirname(library_filename)

    name = config["name"]
    fragments_filename = joinpath(wd, config["fragments"])
    cocktails_filename = joinpath(wd, config["cocktails"])

    fragments = parsefragments(fragments_filename)
    cocktails = parsecocktails(cocktails_filename)

    Library(name, fragments, cocktails)
end


function parsefragments(fragments_filename)
    fragment_data = parsecsv(fragments_filename) # list of (id, smiles, conc, peaks)

    # populate dictionary keyed with cocktail ids
    fragments = Dict{String,LibraryFragment}()
    for fields in fragment_data
        fragment_id = fields[1]
        smiles = fields[2]
        conc = parse(Float64, fields[3])
        shifts = parse.(Float64, fields[4:2:end])
        peakheights = parse.(Float64, fields[5:2:end])
        peak_ids = map(1:length(shifts)) do i
            fragment_id * ('a' + i - 1)
        end
        fragment = LibraryFragment(fragment_id, smiles, conc, peak_ids, shifts, peakheights)
        push!(fragments, fragment_id => fragment)
    end

    return fragments
end


function parsecocktails(filename)
    cocktail_data = parsecsv(filename) # list of (id, name, fragments)

    # populate dictionary keyed with cocktail ids
    cocktails = Dict{String,LibraryCocktail}()
    for fields in cocktail_data
        id = fields[1]
        name = fields[2]
        fragments = fields[3:end]
        cocktail = LibraryCocktail(id, name, fragments)
        push!(cocktails, id => cocktail)
    end

    return cocktails
end


function parsecsv(filename::String)
    result = []
    
    open(filename) do file
        readline(file)      # Skip header line
        
        for line in eachline(file)
            fields = split(line, ',')   # Split line into fields
            push!(result, fields)       # Add to result
        end
    end
    
    return result
end
