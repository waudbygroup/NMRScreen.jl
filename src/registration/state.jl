function initialisestate(config, library, regcocktails)
    state = Dict{String,Any}()
    state["config"] = config
    state["library"] = library
    state["cocktails"] = regcocktails

    state["n_cocktails"] = length(regcocktails)
    state["cocktail_ids"] = [ regcocktail.id for regcocktail in regcocktails ]
    state["current_cocktail_number"] = Observable(1)
    state["current_cocktail"] = lift(i->state["cocktails"][i], state["current_cocktail_number"])

    state["max_peaks"] = maximum([length(cocktail.peak_ids) for cocktail in regcocktails])
    state["n_peaks"] = lift(i -> length(i.peak_ids), state["current_cocktail"])
    state["current_peak_number"] = Observable(1)

    state["reference_spec"] = lift(i -> i.refspec, state["current_cocktail"])
    state["bound_spec"] = lift(i -> i.boundspec, state["current_cocktail"])
    state["reference_dx"] = lift(spec -> abs(data(spec, F1Dim)[2] - data(spec, F1Dim)[1]), state["reference_spec"])
    state["bound_dx"] = lift(spec -> abs(data(spec, F1Dim)[2] - data(spec, F1Dim)[1]), state["bound_spec"])

    # peak positions for current cocktail
    state["current_cocktail_id"] = lift(c -> c.id, state["current_cocktail"])
    state["peak_ids"] = lift(c -> c.peak_ids, state["current_cocktail"])
    state["smiles"] = lift(c -> c.smiles, state["current_cocktail"])
    state["library_x"] = lift(c -> c.library_shifts, state["current_cocktail"])
    state["ref_x"] = @nested_observable state["current_cocktail"] ref_shifts
    state["bound_x"] = @nested_observable state["current_cocktail"] bound_shifts
    state["good"] = @nested_observable_bool state["current_cocktail"] good

    # peak id and smiles for current peak
    state["current_peak_id"] = lift((c, i) -> c[i], state["peak_ids"], state["current_peak_number"])
    state["current_smiles"] = lift((c, i) -> c[i], state["smiles"], state["current_peak_number"])

    # peak heights for current cocktail
    state["ref_y"] = lift(state["ref_x"], state["reference_spec"]) do xs, spec
        [spec[Near(x)]/scale(spec) for x in xs]
    end
    state["bound_y"] = lift(state["bound_x"], state["bound_spec"]) do xs, spec
        [spec[Near(x)]/scale(spec) for x in xs]
    end
    # peak position/heights for current peak
    state["current_library_x"] = lift((c, i) -> c[i], state["library_x"], state["current_peak_number"])
    state["current_ref_x"] = lift((c, i) -> c[i], state["ref_x"], state["current_peak_number"])
    state["current_bound_x"] = lift((c, i) -> c[i], state["bound_x"], state["current_peak_number"])
    state["current_ref_y"] = lift((ys, i) -> ys[i], state["ref_y"], state["current_peak_number"])
    state["current_bound_y"] = lift((ys, i) -> ys[i], state["bound_y"], state["current_peak_number"])
    state["current_good"] = lift((g, i) -> g[i] == 0 ? false : true, state["good"], state["current_peak_number"])

    # points for plotting
    state["ref_points"] = lift(state["ref_y"]) do ys
        # only lift on y because they are updated after x and this will avoid getting out of sync
        [Point2f(x, y) for (x, y) in zip(state["ref_x"][], ys)]
    end
    state["bound_points"] = lift(state["bound_y"]) do ys
        [Point2f(x, y) for (x, y) in zip(state["bound_x"][], ys)]
    end
    state["current_points"] = lift(state["ref_points"], state["bound_points"], state["current_peak_number"]) do ref, bound, i
        [ref[i], bound[i]]
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

    # info
    state["labeltext"] = lift(state["current_cocktail_id"], state["current_peak_id"]) do cocktail, peak
        "Cocktail: $cocktail\nPeak: $peak"
    end

    state["structure"] = lift(i -> smilestoimage(i), state["current_smiles"])

    state["should_close"] = Observable(false)

    state["heatmap_csps1"] = Observable(zeros(state["max_peaks"], state["n_cocktails"]))
    state["heatmap_csps2"] = Observable(zeros(state["max_peaks"], state["n_cocktails"]))
    updateheatmaps!(state)
    state["heatmap_point"] = lift(state["current_peak_number"], state["current_cocktail_number"]) do i, j
        Point2f(i, j)
    end

    return state
end


function updateheatmaps!(state)
    heatmap_csps1 = state["heatmap_csps1"][] # library - reference
    heatmap_csps2 = state["heatmap_csps2"][] # reference - bound
    for i in 1:state["n_cocktails"]
        for j in 1:state["max_peaks"]
            if j > length(state["cocktails"][i].peak_ids)
                heatmap_csps1[j, i] = NaN
                heatmap_csps2[j, i] = NaN
                continue
            end
            libx = state["cocktails"][i].library_shifts[j]
            refx = state["cocktails"][i].ref_shifts[][j]
            boundx = state["cocktails"][i].bound_shifts[][j]
            heatmap_csps1[j, i] = abs(libx - refx)
            heatmap_csps2[j, i] = abs(refx - boundx)
        end
    end
    state["heatmap_csps1"][] = heatmap_csps1
    state["heatmap_csps2"][] = heatmap_csps2
end