function recreatecocktails(state)
    cocktails = map(state["cocktails"]) do cocktail
        peaks = Vector{BasicPeak}()
        for (i, peak_id) in enumerate(cocktail.peak_ids)
            cocktail.good[][i] || continue
            peak = BasicPeak(peak_id,
                cocktail.refspec,
                cocktail.boundspec,
                # cocktail.fragment_ids[i],
                cocktail.id,
                cocktail.smiles[i],
                cocktail.library_shifts[i],
                cocktail.ref_shifts[][i],
                cocktail.bound_shifts[][i])
            push!(peaks, peak)
        end
        Cocktail(cocktail.id, cocktail.name, cocktail.fragment_ids, cocktail.refspec, cocktail.boundspec, peaks)
    end
    return cocktails
end
