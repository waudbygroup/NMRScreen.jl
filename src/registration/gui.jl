function gui!(state)
    fig = Figure(size = (1200, 660))
    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    c3 = Makie.wong_colors()[3]
    c4 = Makie.wong_colors()[4]
    c5 = Makie.wong_colors()[5]
    c6 = Makie.wong_colors()[6]
    c7 = Makie.wong_colors()[7]

    # Top panel
    top_panel = fig[1, 1] = GridLayout()
    button_left = top_panel[1, 1] = Button(fig, label = "<")
    button_right = top_panel[1, 2] = Button(fig, label = ">")
    slider_cocktail_label = top_panel[1,3] = Label(fig, "Cocktail\n(Up/Down)")
    slider_cocktails = top_panel[1,4] = Slider(fig, range = 1:state["n_cocktails"])
    state["slider_cocktails"] = slider_cocktails
    connect!(state["current_cocktail_number"], slider_cocktails.value)
    slider_peak_label = top_panel[1,5] = Label(fig, "Peak\n(Left/Right)")
    slider_peaks = top_panel[1,6] = Slider(fig, range = @lift(1:$(state["n_peaks"])))
    state["slider_peaks"] = slider_peaks
    connect!(state["current_peak_number"], slider_peaks.value)

    cb_label = top_panel[1,7] = Label(fig, "Active\n(SPACE to toggle)")
    cbgood = Checkbox(top_panel[1, 8], checked = state["current_good"][])
    # connect!(cbgood.checked, state["current_good"])
    # when current peak or cocktail changes, update checkbox
    on(state["current_peak_number"]) do _
        cbgood.checked[] = state["current_good"][]
    end
    on(state["current_cocktail_number"]) do _
        cbgood.checked[] = state["current_good"][]
    end
    on(cbgood.checked) do _
        state["good"][][state["current_peak_number"][]] = cbgood.checked[]
        notify(state["good"])
    end

    info_label = top_panel[1,9] = Label(fig, state["labeltext"])
    ax_structure = Axis(top_panel[1,10], aspect = DataAspect(), backgroundcolor = :white, height=150,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,)
    hidedecorations!(ax_structure)
    hidespines!(ax_structure)
    image!(ax_structure, state["structure"])
    button = top_panel[1,11] = Button(fig, label = "Proceed")
    on(button.clicks) do _
        state["should_close"][] = true
    end
    # colsize!(top_panel, 7, Relative(1/3))

    on(button_left.clicks) do _
        left!(state)
    end
    on(button_right.clicks) do _
        right!(state)
    end


    # Middle panel
    # middle_panel = fig[2, 1] = GridLayout()
    spectra_plot = Axis(fig[2,1],
        xreversed=true,
        yzoomlock=true,
        ypanlock=true,
        xrectzoom=true,
        yrectzoom=false,
        xlabel="Chemical shift (ppm)",
        ylabel="Intensity",
        )
    hideydecorations!(spectra_plot)

    state["current_color"] = lift(good -> good ? c3 : c6, state["current_good"])

    vlines!(spectra_plot, state["library_x"], color=:grey, ymin=0.9, label="Library")
    vlines!(spectra_plot, state["current_library_x"], color=state["current_color"], ymin=0.85, linewidth=3)
    lines!(spectra_plot, state["reference_plot"], label = "Reference", color=c5)
    lines!(spectra_plot, state["bound_plot"], label = "Bound", color=c7)
    axislegend(spectra_plot)
    
    scatter!(spectra_plot, state["ref_points"], color=c1, markersize=12)
    scatter!(spectra_plot, state["bound_points"], color=c2, markersize=12)
    scatter!(spectra_plot, state["current_points"], markersize=15, color=state["current_color"])


    # Bottom panel
    bottom_panel = fig[3,1] = GridLayout()

    ax_heatmap1 = bottom_panel[1,1] = Axis(fig,
        xlabel="Peak",
        ylabel="Cocktail",
        yreversed=true,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,
        backgroundcolor=:grey10,
        title="Reference vs library peak positions",
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm1 = heatmap!(ax_heatmap1, state["heatmap_csps1"],
        colorrange=lift(z->(0, maximum(filter(!isnan, z))), state["heatmap_csps1"]))
    cb = Colorbar(bottom_panel[1,2], hm1, label="Chemical shift difference (ppm)")
    scatter!(ax_heatmap1, state["heatmap_point"], color=c4, marker=:star5, markersize=12)

    ax_heatmap2 = bottom_panel[1,3] = Axis(fig,
        xlabel="Peak",
        ylabel="Cocktail",
        yreversed=true,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,
        backgroundcolor=:grey10,
        title="Bound vs reference peak positions",
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
        
    hm2 = heatmap!(ax_heatmap2, state["heatmap_csps2"],
        colorrange=lift(z->(0, maximum(filter(!isnan, z))), state["heatmap_csps2"]))
    cb = Colorbar(bottom_panel[1,4], hm2, label="Chemical shift difference (ppm)")
    scatter!(ax_heatmap2, state["heatmap_point"], color=c4, marker=:star5, markersize=12)
    # Event handler for heatmap click
    on(events(ax_heatmap1).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            plt, i = pick(fig)
            plt == hm1 || return
            i > 0 || return
            idx = CartesianIndices(state["heatmap_csps1"][])[i]
            peak_idx = idx[1]
            cocktail_idx = idx[2]
            if peak_idx <= length(state["cocktails"][cocktail_idx].peak_ids)
                set_close_to!(slider_cocktails, cocktail_idx)
                set_close_to!(slider_peaks, peak_idx)
            end
        end
    end
    on(events(ax_heatmap2).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            plt, i = pick(fig)
            plt == hm2 || return
            i > 0 || return
            idx = CartesianIndices(state["heatmap_csps2"][])[i]
            peak_idx = idx[1]
            cocktail_idx = idx[2]
            if peak_idx <= length(state["cocktails"][cocktail_idx].peak_ids)
                set_close_to!(slider_cocktails, cocktail_idx)
                set_close_to!(slider_peaks, peak_idx)
            end
        end
    end


    # on(slider_peaks.value) do n
    on(state["current_ref_x"]) do x
        xl1, xl2 = spectra_plot.xaxis.attributes.limits[]
        mn,mx = extrema(data(state["reference_spec"][],F1Dim))
        sw = mx - mn
        xw = xl2 - xl1
        # x = state["ref_points"][][n][1]
        # x = p.reference_shift
        nx1 = x - xw / 2
        nx2 = x + xw / 2
        if sw > xw # zoomed in
            if nx1 < mn
                nx1 = mn
                nx2 = mn + xw
            elseif nx2 > mx
                nx2 = mx
                nx1 = mx - xw
            end
        else # zoomed out
            if nx1 > mn
                nx1 = mn
                nx2 = mn + xw
            elseif nx2 < mx
                nx2 = mx
                nx1 = mx - xw
            end
        end
        xlims!(spectra_plot, (nx2, nx1))
    end
    on(state["reference_plot"]) do spec
        # reset y limits
        mn,mx = extrema([p[2] for p in spec])
        ylims!(spectra_plot, (-0.1mx, 1.1mx))
    end
    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if ispressed(fig, Exclusively(Keyboard.left))
                left!(state)
            elseif ispressed(fig, Exclusively(Keyboard.right))
                right!(state)
            elseif ispressed(fig, Keyboard.left & (Keyboard.left_shift | Keyboard.right_shift))
                navigate_peaks_within_fragment!(state, -1)
            elseif ispressed(fig, Keyboard.right & (Keyboard.left_shift | Keyboard.right_shift))
                navigate_peaks_within_fragment!(state, 1)
            elseif event.key == Keyboard.down
                i = state["current_cocktail_number"][]
                if i < state["n_cocktails"]
                    set_close_to!(slider_cocktails, i+1)
                end
            elseif event.key == Keyboard.up
                i = state["current_cocktail_number"][]
                if i > 1
                    set_close_to!(slider_cocktails, i-1)
                end
            elseif ispressed(fig, Exclusively(Keyboard.left_bracket))
                nudge_ref_left!(state)
            elseif ispressed(fig, Exclusively(Keyboard.right_bracket))
                nudge_ref_right!(state)
            elseif ispressed(fig, Exclusively(Keyboard.comma))
                nudge_bound_left!(state)
            elseif ispressed(fig, Exclusively(Keyboard.period))
                nudge_bound_right!(state)
            elseif ispressed(fig, Keyboard.left_bracket & (Keyboard.left_shift | Keyboard.right_shift))
                # Find next peak in reference spectrum
                current_pos = state["ref_x"][][state["current_peak_number"][]]
                new_pos = find_nearest_peak(state["reference_spec"][], current_pos, -1)
                state["ref_x"][][state["current_peak_number"][]] = new_pos
                notify(state["ref_x"])
                updateheatmaps!(state)
            elseif ispressed(fig, Keyboard.right_bracket & (Keyboard.left_shift | Keyboard.right_shift))
                # Find next peak in reference spectrum
                current_pos = state["ref_x"][][state["current_peak_number"][]]
                new_pos = find_nearest_peak(state["reference_spec"][], current_pos, +1)
                state["ref_x"][][state["current_peak_number"][]] = new_pos
                notify(state["ref_x"])
                updateheatmaps!(state)
            elseif ispressed(fig, Keyboard.comma & (Keyboard.left_shift | Keyboard.right_shift))
                # Find next peak in bound spectrum
                current_pos = state["bound_x"][][state["current_peak_number"][]]
                new_pos = find_nearest_peak(state["bound_spec"][], current_pos, -1)
                state["bound_x"][][state["current_peak_number"][]] = new_pos
                notify(state["bound_x"])
                updateheatmaps!(state)
            elseif ispressed(fig, Keyboard.period & (Keyboard.left_shift | Keyboard.right_shift))
                # Find next peak in bound spectrum
                current_pos = state["bound_x"][][state["current_peak_number"][]]
                new_pos = find_nearest_peak(state["bound_spec"][], current_pos, +1)
                state["bound_x"][][state["current_peak_number"][]] = new_pos
                notify(state["bound_x"])
                updateheatmaps!(state)
            end
            if event.action == Keyboard.press && event.key == Keyboard.space
                cbgood.checked[] = !cbgood.checked[]
            end
        end
    end

    rowsize!(fig.layout, 1, Relative(0.2))
    rowsize!(fig.layout, 3, Relative(0.4))

    showhelp()

    display(fig, float = false)
    while !state["should_close"][]
        sleep(0.1)
        if !isopen(fig.scene)
            break
        end
    end

    state["should_close"][] = false # reset for next iteration
    GLMakie.closeall()
end

    
function showhelp()
    @info """
# NMR Cocktail Peak Alignment

Keyboard shortcuts:
- ←/→: Navigate between peaks
- SHIFT + ←/→: Navigate between peaks within the same fragment
- ↑/↓: Navigate between cocktails
- '[' / ']': Nudge reference peak left/right (SHIFT for next maximum)
- ',' / '.': Nudge bound peak left/right (SHIFT for next maximum)
- Control-click: Reset zoom
- Right-drag: Pan spectrum
    """
end

function nudge_ref_left!(state, n=1)
    state["ref_x"][][state["current_peak_number"][]] += n * state["reference_dx"][]
    notify(state["ref_x"])
    updateheatmaps!(state)
end
function nudge_ref_right!(state, n=1)
    state["ref_x"][][state["current_peak_number"][]] -= n * state["reference_dx"][]
    notify(state["ref_x"])
    updateheatmaps!(state)
end
function nudge_bound_left!(state, n=1)
    state["bound_x"][][state["current_peak_number"][]] += n * state["bound_dx"][]
    notify(state["bound_x"])
    updateheatmaps!(state)
end
function nudge_bound_right!(state, n=1)
    state["bound_x"][][state["current_peak_number"][]] -= n * state["bound_dx"][]
    notify(state["bound_x"])
    updateheatmaps!(state)
end

function left!(state)
    i = state["current_cocktail_number"][]
    j = state["current_peak_number"][]
    if j > 1
        j -= 1
        set_close_to!(state["slider_peaks"], j)
    elseif i > 1
        i -= 1
        j = length(state["cocktails"][i].peak_ids)
        set_close_to!(state["slider_cocktails"], i)
        set_close_to!(state["slider_peaks"], j)
    end
end
function right!(state)
    i = state["current_cocktail_number"][]
    j = state["current_peak_number"][]
    if j < length(state["cocktails"][i].peak_ids)
        j += 1
        set_close_to!(state["slider_peaks"], j)
    elseif i < state["n_cocktails"]
        i += 1
        j = 1
        set_close_to!(state["slider_cocktails"], i)
        set_close_to!(state["slider_peaks"], j)
    end
end

function navigate_peaks_within_fragment!(state, direction)
    current_cocktail = state["current_cocktail"][]
    current_peak_number = state["current_peak_number"][]
    current_fragment_id = current_cocktail.fragment_ids[current_peak_number]

    # find other peaks matching this fragment id, then set the current peak number to the previous/next one (depending on direction)
    other_peak_numbers = findall(==(current_fragment_id), current_cocktail.fragment_ids)

    # Find the current index in other_peak_numbers
    current_index = findfirst(==(current_peak_number), other_peak_numbers)
    
    if isnothing(current_index)
        error("Current peak number not found in fragment peaks")
    end
    
    # Calculate the new index based on direction
    if direction == 1
        new_index = current_index % length(other_peak_numbers) + 1
    elseif direction == -1
        new_index = (current_index - 2 + length(other_peak_numbers)) % length(other_peak_numbers) + 1
    else
        error("Invalid direction: $direction. Must be +/- 1")
    end
    
    # Update the current peak number in the state
    set_close_to!(state["slider_peaks"], other_peak_numbers[new_index])
end

function find_nearest_peak(spec, current_pos, direction; snr_threshold=3)
    x = data(spec, F1Dim)
    y = data(spec[:,1]) / spec[:noise]

    pks = findmaxima(y)     # detect peaks
    peakheights!(pks; min=snr_threshold) # filter list

    peak_shifts = x[pks.indices]
    
    # Find peaks in the specified direction from current position
    if direction > 0
        candidates = peak_shifts[peak_shifts .< current_pos]
        isempty(candidates) && return current_pos
        return maximum(candidates)
    else
        candidates = peak_shifts[peak_shifts .> current_pos]
        isempty(candidates) && return current_pos
        return minimum(candidates)
    end
end