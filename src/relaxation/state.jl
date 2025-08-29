function initialisestate(config, library, cocktails)
    state = Dict{String,Any}()
    state["config"] = config
    state["library"] = library
    state["cocktails"] = cocktails

    state["n_cocktails"] = length(cocktails)
    state["cocktail_ids"] = [ cocktail.id for cocktail in cocktails ]
    state["current_cocktail_number"] = Observable(1)
    state["current_cocktail"] = lift(i->state["cocktails"][i], state["current_cocktail_number"])

    state["max_peaks"] = maximum([length(cocktail.peaks) for cocktail in cocktails])
    state["n_peaks"] = lift(i -> length(i.peaks), state["current_cocktail"])
    state["current_peak_number"] = Observable(1)
    state["current_peak"] = lift((i, j) -> i.peaks[j], state["current_cocktail"], state["current_peak_number"])

    state["library_shifts"] = lift(i -> [p.library_shift for p in i.peaks], state["current_cocktail"])
    state["ref_points"] = lift(i -> [Point2f(p.reference_shift, p.reference_intensity) for p in i.peaks], state["current_cocktail"])
    state["bound_points"] = lift(i -> [Point2f(p.bound_shift, p.bound_intensity) for p in i.peaks], state["current_cocktail"])
    state["current_library_x"] = lift(p -> p.library_shift, state["current_peak"])
    state["current_ref_x"] = lift(p -> p.reference_shift, state["current_peak"])
    state["current_ref_y"] = lift(p -> p.reference_intensity, state["current_peak"])
    state["current_bound_x"] = lift(p -> p.bound_shift, state["current_peak"])
    state["current_bound_y"] = lift(p -> p.bound_intensity, state["current_peak"])

    state["reference_spec"] = lift(i -> i.refspec, state["current_cocktail"])
    state["bound_spec"] = lift(i -> i.boundspec, state["current_cocktail"])

    state["reference_relaxation"] = lift(i -> i.reference_relaxation, state["current_peak"])
    state["bound_relaxation"] = lift(i -> i.bound_relaxation, state["current_peak"])

    state["refR2"] = lift(i -> refR2(i), state["current_peak"])
    state["DeltaR2"] = lift(i -> DeltaR2(i), state["current_peak"])
    state["II0"] = lift(i -> II0(i), state["current_peak"])
    state["csp"] = lift(i -> csp(i), state["current_peak"])

    state["csps_heatmap"] = Observable(csps(state))
    state["II0_heatmap"] = Observable(II0s(state))
    state["DeltaR2_heatmap"] = Observable(DeltaR2s(state))
    state["refR2_heatmap"] = Observable(refR2s(state))
    state["reducedchi2_heatmap"] = Observable(reducedchi2s(state))

    state["umap_xs"] = Observable(umap_xs(state))
    state["umap_ys"] = Observable(umap_ys(state))
    state["umap_csps"] = Observable(umap_csps(state))
    state["umap_II0s"] = Observable(umap_II0s(state))
    state["umap_DeltaR2s"] = Observable(umap_DeltaR2s(state))
    state["umap_refR2s"] = Observable(umap_refR2s(state))
    state["umap_reducedchi2s"] = Observable(umap_reducedchi2s(state))

    state["smiles"] = lift(i -> library.fragments[i.fragment_id].smiles, state["current_peak"])
    state["structure"] = lift(i -> smilestoimage(i), state["smiles"])

    addguistuff!(state)
    return state
end

function addguistuff!(state)
    state["nfragmentpeaks"] = lift(state["current_cocktail"], state["current_peak_number"]) do current_cocktail, current_peak_number
        current_fragment_id = current_cocktail.fragment_ids[current_peak_number]
        other_peaks = findall(==(current_fragment_id), current_cocktail.fragment_ids)
        current_index = findfirst(==(current_peak_number), other_peaks)
        (current_index, length(other_peaks))
    end
    state["labeltext"] = lift(state["current_peak"], state["nfragmentpeaks"]) do peak, nfragmentpeaks
        current, n = nfragmentpeaks
        if n == 1
            "Cocktail: $(peak.cocktail_id)\nPeak: $(peak.peak_id)"
        else
            "Cocktail: $(peak.cocktail_id)\nPeak: $(peak.peak_id)\n($current of $n)"
        end
    end
    # state["labeltext"] = lift(state["current_peak"]) do peak
    #     "Cocktail: $(peak.cocktail_id)\nPeak: $(peak.peak_id)"
    # end
    state["peakinfotext"] = lift(state["current_peak"]) do peak
        # Cocktail: $(peak.cocktail_id), Peak: $(peak.peak_id)
        """I/I₀: $(II0(peak)), Δδ: $(round(csp(peak), digits=2)) ppm
        Ref. R₂: $(refR2(peak)) s-1, ΔR₂: $(DeltaR2(peak)) s-1
        """
    end
    state["II0_csp_text"] = lift(state["current_peak"]) do peak
        "I/I₀: $(II0(peak)), Δδ: $(round(csp(peak), digits=2)) ppm"
    end
    state["R2_text"] = lift(state["current_peak"]) do peak
        "Ref. R₂: $(refR2(peak)) s-1, ΔR₂: $(DeltaR2(peak)) s-1"
    end

    # plotting data
    state["reference_plot"] = lift(state["reference_spec"]) do spec
        x = data(spec, F1Dim)
        y = data(spec[:,1]) / scale(spec)
        [Point2f(p...) for p in zip(x, y)]
    end
    state["bound_plot"] = lift(state["bound_spec"]) do spec
        x = data(spec, F1Dim)
        y = data(spec[:,1]) / scale(spec)
        [Point2f(p...) for p in zip(x, y)]
    end
    
    state["heatmap_menu"] = Observable("DeltaR2")
    state["heatmap_data"] = Observable(randn(state["max_peaks"], state["n_cocktails"]))
    state["heatmap_label"] = Observable("ΔR₂ (s⁻¹)")
    state["heatmap_limits"] = Observable{Any}((0, 20.))
    state["heatmap_cm"] = Observable{Any}(:viridis)
    state["heatmap_point"] = lift(state["current_peak_number"], state["current_cocktail_number"]) do i, j
        Point2f(i, j)
    end

    state["umap_color"] = Observable(randn(sum([length(cocktail.peaks) for cocktail in state["cocktails"]])))
    state["umap_current_x"] = lift(i -> i.umap_x, state["current_peak"])
    state["umap_current_y"] = lift(i -> i.umap_y, state["current_peak"])

    state["should_close"] = Observable(false)
end

function refR2s(state)
    refR2s = fill(NaN, state["max_peaks"], state["n_cocktails"])
    for (i, cocktail) in enumerate(state["cocktails"])
        for (j, peak) in enumerate(cocktail.peaks)
            refR2s[j, i] = Measurements.value(refR2(peak))
        end
    end

    return refR2s
end

function DeltaR2s(state)
    deltaR2s = fill(NaN, state["max_peaks"], state["n_cocktails"])
    for (i, cocktail) in enumerate(state["cocktails"])
        for (j, peak) in enumerate(cocktail.peaks)
            deltaR2s[j, i] = Measurements.value(DeltaR2(peak))
        end
    end

    return deltaR2s
end

function II0s(state)
    II0s = fill(NaN, state["max_peaks"], state["n_cocktails"])
    for (i, cocktail) in enumerate(state["cocktails"])
        for (j, peak) in enumerate(cocktail.peaks)
            II0s[j, i] = Measurements.value(II0(peak))
        end
    end

    return II0s
end

function csps(state)
    csps = fill(NaN, state["max_peaks"], state["n_cocktails"])
    for (i, cocktail) in enumerate(state["cocktails"])
        for (j, peak) in enumerate(cocktail.peaks)
            csps[j, i] = abs(Measurements.value(csp(peak)))
        end
    end

    return csps
end

function reducedchi2s(state)
    reducedchi2s = fill(NaN, state["max_peaks"], state["n_cocktails"])
    for (i, cocktail) in enumerate(state["cocktails"])
        for (j, peak) in enumerate(cocktail.peaks)
            reducedchi2s[j, i] = reducedchi2(peak)
        end
    end

    return reducedchi2s
end

function umap_xs(state)
    umap_xs = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_xs, peak.umap_x)
        end
    end
    return umap_xs
end

function umap_ys(state)
    umap_ys = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_ys, peak.umap_y)
        end
    end
    return umap_ys
end

function umap_refR2s(state)
    umap_refR2s = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_refR2s, Measurements.value(refR2(peak)))
        end
    end
    return umap_refR2s
end

function umap_DeltaR2s(state)
    umap_DeltaR2s = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_DeltaR2s, Measurements.value(DeltaR2(peak)))
        end
    end
    return umap_DeltaR2s
end

function umap_II0s(state)
    umap_II0s = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_II0s, Measurements.value(II0(peak)))
        end
    end
    return umap_II0s
end

function umap_csps(state)
    umap_csps = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_csps, abs(Measurements.value(csp(peak))))
        end
    end
    return umap_csps
end

function umap_reducedchi2s(state)
    umap_reducedchi2s = []
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            push!(umap_reducedchi2s, reducedchi2(peak))
        end
    end
    return umap_reducedchi2s
end