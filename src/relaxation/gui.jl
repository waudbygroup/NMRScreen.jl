function gui!(state)
    fig = Figure(size = (1200, 800))
    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    c3 = Makie.wong_colors()[3]
    c4 = Makie.wong_colors()[6]

    # Top panel
    top_panel = fig[1, 1] = GridLayout()
    button_left = top_panel[1, 1] = Button(fig, label = "<")
    button_right = top_panel[1, 2] = Button(fig, label = ">")
    slider_cocktail_label = top_panel[1,3] = Label(fig, "Cocktail")
    slider_cocktails = top_panel[1,4] = Slider(fig, range = 1:state["n_cocktails"])
    connect!(state["current_cocktail_number"], slider_cocktails.value)
    slider_peak_label = top_panel[1,5] = Label(fig, "Peak")
    slider_peaks = top_panel[1,6] = Slider(fig, range = @lift(1:$(state["n_peaks"])))
    state["slider_peaks"] = slider_peaks
    connect!(state["current_peak_number"], slider_peaks.value)
    info_label = top_panel[1,7] = Label(fig, state["labeltext"])
    ax_structure = Axis(top_panel[1,8], aspect = DataAspect(), backgroundcolor = :white, height=120,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,)
    hidedecorations!(ax_structure)
    hidespines!(ax_structure)
    image!(ax_structure, state["structure"])
    back_button = top_panel[1,9] = Button(fig, label = "Back")
    on(back_button.clicks) do _
        go_back(state)
    end
    save_button = top_panel[1,10] = Button(fig, label = "Save")
    on(save_button.clicks) do _
        write_results(state)
    end
    quit_button = top_panel[1,11] = Button(fig, label = "Quit")
    on(quit_button.clicks) do _
        quit(state)
    end
    colsize!(top_panel, 7, Relative(1/3))

    on(button_left.clicks) do _
        i = state["current_cocktail_number"][]
        j = state["current_peak_number"][]
        if j > 1
            j -= 1
            set_close_to!(slider_peaks, j)
        elseif i > 1
            i -= 1
            j = length(state["cocktails"][i].peaks)
            set_close_to!(slider_cocktails, i)
            set_close_to!(slider_peaks, j)
        end
    end
    on(button_right.clicks) do _
        i = state["current_cocktail_number"][]
        j = state["current_peak_number"][]
        if j < length(state["cocktails"][i].peaks)
            j += 1
            set_close_to!(slider_peaks, j)
        elseif i < state["n_cocktails"]
            i += 1
            j = 1
            set_close_to!(slider_cocktails, i)
            set_close_to!(slider_peaks, j)
        end
    end


    # Middle panel
    # middle_panel = fig[2, 1] = GridLayout()
    ref_spectra_plot = Axis(fig[2,1],
        xreversed=true,
        yzoomlock=true,
        ypanlock=true,
        xrectzoom=true,
        yrectzoom=false,
        xlabel="Chemical shift (ppm)",
        ylabel="Intensity",
        )
    hideydecorations!(ref_spectra_plot)

    vlines!(ref_spectra_plot, state["library_shifts"], color=c3, ymin=0.9, label="Library")
    vlines!(ref_spectra_plot, state["current_library_x"], color=c4, ymin=0.85, linewidth=3)
    lines!(ref_spectra_plot, state["reference_plot"], label = "Reference", color=c1)
    lines!(ref_spectra_plot, state["bound_plot"], label = "Bound", color=c2)
    axislegend(ref_spectra_plot)
    
    scatter!(ref_spectra_plot, state["ref_points"], color=c1)
    scatter!(ref_spectra_plot, state["bound_points"], color=c2)
    scatter!(ref_spectra_plot, state["current_ref_x"], state["current_ref_y"], color = c4, markersize=12)
    scatter!(ref_spectra_plot, state["current_bound_x"], state["current_bound_y"], color = c4, markersize=12)

    # on(slider_peaks.value) do n
    on(state["current_peak"]) do p
        xl1, xl2 = ref_spectra_plot.xaxis.attributes.limits[]
        mn,mx = extrema(data(state["reference_spec"][],F1Dim))
        sw = mx - mn
        xw = xl2 - xl1
        # x = state["ref_points"][][n][1]
        x = p.reference_shift
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
        xlims!(ref_spectra_plot, (nx2, nx1))
    end
    on(state["reference_plot"]) do spec
        # reset y limits
        mn,mx = extrema([p[2] for p in spec])
        ylims!(ref_spectra_plot, (-0.1mx, 1.1mx))
    end


    # Bottom-left panel
    bottom_panel = fig[3, 1] = GridLayout()
    peak_info_label = bottom_panel[1,1] = Label(fig, state["peakinfotext"])
    
    # setup relaxation plot
    relaxation_plot = bottom_panel[2,1] = Axis(fig,
        xlabel="Relaxation time (ms)",
        ylabel="Intensity",
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,
        limits=(nothing, (-.1, 1.1)))
    lines!(relaxation_plot, @lift($(state["reference_relaxation"]).typred), color=c1)
    lines!(relaxation_plot, @lift($(state["bound_relaxation"]).typred), color=c2)
    errorbars!(relaxation_plot, @lift($(state["reference_relaxation"]).tye), color=c1)
    errorbars!(relaxation_plot, @lift($(state["bound_relaxation"]).tye), color=c2)
    scatter!(relaxation_plot, @lift($(state["reference_relaxation"]).ty), label = "Reference", color=c1)
    scatter!(relaxation_plot, @lift($(state["bound_relaxation"]).ty), label = "Bound", color=c2)
    axislegend(relaxation_plot)
    

    # # Bottom-right panel
    # data_menu = menus[1,1] = Menu(fig, options = ["DeltaR2", "Relative Intensity", "CSP", "Base R2"])
    menus = bottom_panel[1,2] = GridLayout()
    data_menu = menus[1,1] = Menu(fig, options = ["ΔR₂", "Relative Intensity", "Chemical Shift Perturbation", "Reference R₂", "Reduced χ²"])
    on(data_menu.selection) do s
        if s == "ΔR₂"
            connect!(state["heatmap_data"], state["DeltaR2_heatmap"])
            connect!(state["umap_color"], state["umap_DeltaR2s"])
            state["heatmap_label"][] = "ΔR₂ (s⁻¹)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["DeltaR2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Relative Intensity"
            connect!(state["heatmap_data"], state["II0_heatmap"])
            connect!(state["umap_color"], state["umap_II0s"])
            state["heatmap_label"][] = "I/I₀"
            state["heatmap_limits"][] = (0., 1.)
            state["heatmap_cm"][] = Reverse(:viridis)
        elseif s == "Reference R₂"
            connect!(state["heatmap_data"], state["refR2_heatmap"])
            connect!(state["umap_color"], state["umap_refR2s"])
            state["heatmap_label"][] = "Reference R₂ (s⁻¹)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["refR2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Chemical Shift Perturbation"
            connect!(state["heatmap_data"], state["csps_heatmap"])
            connect!(state["umap_color"], state["umap_csps"])
            state["heatmap_label"][] = "Chemical Shift Perturbation (ppm)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["csps_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Reduced χ²"
            connect!(state["heatmap_data"], state["reducedchi2_heatmap"])
            connect!(state["umap_color"], state["umap_reducedchi2s"])
            state["heatmap_label"][] = "Reduced χ²"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["reducedchi2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        end
    end
    notify(data_menu.selection)
    plot_menu = menus[1,2] = Menu(fig, options = ["Heatmap", "UMAP"])
    # Offscreen box for storing idle plot
    idle_plot = GridLayout(bbox = BBox(-200, -100, 0, 100))

    heatmap = GridLayout()
    ax_heatmap = heatmap[1,1] = Axis(fig,
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
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    heatmap_graph = heatmap!(ax_heatmap, state["heatmap_data"],
        colormap=state["heatmap_cm"],
        colorrange=state["heatmap_limits"])
    scatter!(ax_heatmap, state["heatmap_point"], color=c4, marker=:star5, markersize=12)
    # Event handler for heatmap click
    on(events(ax_heatmap).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            plt, i = pick(fig)
            plt == heatmap_graph || return
            i > 0 || return
            idx = CartesianIndices(state["heatmap_data"][])[i]
            peak_idx = idx[1]
            cocktail_idx = idx[2]
            if peak_idx <= length(state["cocktails"][cocktail_idx].peaks)
                set_close_to!(slider_cocktails, cocktail_idx)
                set_close_to!(slider_peaks, peak_idx)
                # state["current_cocktail_number"][] = cocktail_idx
                # state["current_peak_number"][] = peak_idx
            end
        end
    end
    on(events(ax_heatmap).scroll) do (_, dy)
        s = data_menu.selection[]
        if s == "ΔR₂"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1..100))
        elseif s == "Relative Intensity"
            state["heatmap_limits"][] = (clamp(state["heatmap_limits"][][1] + dy/10, 0..0.9), 1.)
        elseif s == "Reference R₂"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1..100))
        elseif s == "Chemical Shift Perturbation"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy/100, 0.01,2))
        elseif s == "Reduced χ²"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1,100))
        end
    end

    chemical_space = GridLayout()
    ax_umap = chemical_space[1,1] = Axis(fig,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=true,
        yrectzoom=true,
        xpanlock=false,
        ypanlock=false)
    umap_plot = scatter!(ax_umap, 
        state["umap_xs"], 
        state["umap_ys"],
        color=state["umap_color"], 
        colormap=state["heatmap_cm"],
        colorrange=state["heatmap_limits"],
        markersize=12)
    #scatter!(ax_umap, state["umap_current_x"], state["umap_current_y"], color=c4, marker='×', markersize=20)
    on(events(ax_umap).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            plt, i = pick(fig)
            plt == umap_plot || return
            i > 0 || return
            peak_idx = i
            cocktail_idx = 1
            while peak_idx > length(state["cocktails"][cocktail_idx].peaks)
                peak_idx -= length(state["cocktails"][cocktail_idx].peaks)
                cocktail_idx += 1
            end

            if peak_idx <= length(state["cocktails"][cocktail_idx].peaks)
                set_close_to!(slider_cocktails, cocktail_idx)
                set_close_to!(slider_peaks, peak_idx)
                # state["current_cocktail_number"][] = cocktail_idx
                # state["current_peak_number"][] = peak_idx
            end
        end
    end
    on(events(ax_umap).scroll) do (_, dy)
        s = data_menu.selection[]
        if s == "ΔR₂"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1..100))
        elseif s == "Relative Intensity"
            state["heatmap_limits"][] = (clamp(state["heatmap_limits"][][1] + dy/10, 0..0.9), 1.)
        elseif s == "Reference R₂"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1..100))
        elseif s == "Chemical Shift Perturbation"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy/100, 0.01,2))
        elseif s == "Reduced χ²"
            state["heatmap_limits"][] = (0., clamp(state["heatmap_limits"][][2] + dy, 1,100))
        end
    end
    
    bottom_panel[2,2] = heatmap
    on(plot_menu.selection) do s
        if s == "Heatmap"
            bottom_panel[2,2] = heatmap
            idle_plot[1,1] = chemical_space
            heatmap_cb = Colorbar(bottom_panel[2,3], heatmap_graph, label=state["heatmap_label"])
        elseif s == "UMAP"
            bottom_panel[2,2] = chemical_space
            idle_plot[1,1] = heatmap
            umap_cb = Colorbar(bottom_panel[2,3], umap_plot, label=state["heatmap_label"])
        end
    end
    notify(plot_menu.selection)

    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if ispressed(fig, Exclusively(Keyboard.left))
                i = state["current_cocktail_number"][]
                j = state["current_peak_number"][]
                if j > 1
                    j -= 1
                    set_close_to!(slider_peaks, j)
                elseif i > 1
                    i -= 1
                    j = length(state["cocktails"][i].peaks)
                    set_close_to!(slider_cocktails, i)
                    set_close_to!(slider_peaks, j)
                end
            elseif ispressed(fig, Exclusively(Keyboard.right))
                i = state["current_cocktail_number"][]
                j = state["current_peak_number"][]
                if j < length(state["cocktails"][i].peaks)
                    j += 1
                    set_close_to!(slider_peaks, j)
                elseif i < state["n_cocktails"]
                    i += 1
                    j = 1
                    set_close_to!(slider_cocktails, i)
                    set_close_to!(slider_peaks, j)
                end
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
            end
        end
    end

    colsize!(bottom_panel, 1, Relative(1/4))
    rowsize!(fig.layout, 1, Relative(0.12))
    rowsize!(fig.layout, 3, Relative(0.58))

    showhelp()

    state["finished"] = false
    display(fig, float = false)
    while !state["should_close"][]
        sleep(0.1)
        if !isopen(fig.scene)
            break
        end
    end
    
    GLMakie.closeall()

    return state["finished"]
end


function showhelp()
    @info """
# NMR Fragment Screening Analysis

Keyboard shortcuts:
- ←/→: Navigate between peaks
- SHIFT + ←/→: Navigate between peaks within the same fragment
- ↑/↓: Navigate between cocktails
- Mouse wheel: Adjust heatmap scale / spectrum zoom level
- Click points on heatmap / chemical space map: Select peak
- Control-click: Reset zoom of spectrum / chemical space map
- Right-drag: Pan spectrum / chemical space map
    """
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

function go_back(state)
    state["finished"] = false
    state["should_close"][] = true
end

function quit(state)
    state["finished"] = true
    state["should_close"][] = true
end