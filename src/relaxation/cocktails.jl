function prepcocktails(cocktails)
    # input: cocktails where peaks are BasicPeaks
    # output: cocktails where peaks are ScreeningPeaks, also initialize the UMAP coordinates fields
    for cocktail in cocktails
        for (i, peak) in enumerate(cocktail.peaks)
            newpeak = ScreeningPeak(peak)
            cocktail.peaks[i] = newpeak
        end
        sort!(cocktail.peaks, by=sp -> sp.reference_shift, rev=true)
    end
    
    embeddings = get_embeddings(cocktails)
    peak_count = 0
    for cocktail in cocktails
        for (i, peak) in enumerate(cocktail.peaks)
            peak_count += 1
            newpeak = ScreeningPeak(peak, embeddings[1, peak_count], embeddings[2, peak_count])
            cocktail.peaks[i] = newpeak
        end
    end
    return cocktails
end

function get_embeddings(cocktails)
    smiles = []
    DeltaR2s = []
    for cocktail in cocktails
        for peak in cocktail.peaks
            push!(smiles, peak.smiles)
            push!(DeltaR2s, Measurements.value(DeltaR2(peak)))
        end
    end
    fingerprint_matrix = []
    for (i, entry) in enumerate(smiles)
        mol = get_mol(entry)
        fingerprint = get_morgan_fp(mol)
        fingerprint_split = split(fingerprint, "")
        fingerprint_parsed = []
        for bits in fingerprint_split
            push!(fingerprint_parsed, parse(Float64, bits))
        end
        max_DeltaR2 = maximum(DeltaR2s)
        if DeltaR2s[i] < 0.0
            DeltaR2_factor = 0.0
        else
            DeltaR2_factor = (DeltaR2s[i] / max_DeltaR2) * 20.0
        end
        push!(fingerprint_parsed, DeltaR2_factor)
        push!(fingerprint_matrix, fingerprint_parsed)
    end
    fingerprint_matrix = hcat(fingerprint_matrix...)
    return umap(fingerprint_matrix, 2; n_neighbors=20, min_dist=0.1, metric=Jaccard())
end
