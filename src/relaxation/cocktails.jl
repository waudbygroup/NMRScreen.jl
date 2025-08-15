function prepcocktails(cocktails)
    embeddings = get_embeddings(cocktails)
    # input: cocktails where peaks are BasicPeaks
    # output: cocktails where peaks are ScreeningPeaks, also contains UMAP coordinates
    peak_count = 0
    for cocktail in cocktails
        for (i, peak) in enumerate(cocktail.peaks)
            peak_count += 1
            newpeak = ScreeningPeak(peak, embeddings[1, peak_count], embeddings[2, peak_count])
            cocktail.peaks[i] = newpeak
        end
        sort!(cocktail.peaks, by=sp -> sp.reference_shift, rev=true)
    end

    return cocktails
end

function get_embeddings(cocktails)
    smiles = []
    for cocktail in cocktails
        for peak in cocktail.peaks
            push!(smiles, peak.smiles)
        end
    end
    fingerprint_matrix = []
    for i in smiles
        mol = get_mol(i)
        fingerprint = get_morgan_fp(mol)
        fingerprint_split = split(fingerprint, "")
        fingerprint_parsed = []
        for bits in fingerprint_split
            push!(fingerprint_parsed, parse.(Float64, bits))
        end
        push!(fingerprint_matrix, fingerprint_parsed)
    end
    fingerprint_matrix = hcat(fingerprint_matrix...)
    return umap(fingerprint_matrix, 2; n_neighbors=20, min_dist=0.1, metric=Jaccard())
end
