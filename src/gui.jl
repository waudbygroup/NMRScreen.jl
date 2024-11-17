function gui(state)
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
    save_button = top_panel[1,9] = Button(fig, label = "Save")
    on(save_button.clicks) do _
        write_results(state)
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

    vlines!(ref_spectra_plot, state["library_shifts"], color=c3, ymin=0.85, label="Library")
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
    # menu = bottom_panel[1,2] = Menu(fig, options = ["DeltaR2", "Relative Intensity", "CSP", "Base R2"])
    menu = bottom_panel[1,2] = Menu(fig, options = ["ΔR₂", "Relative Intensity", "Chemical Shift Perturbation", "Reference R₂", "Reduced χ²"])
    on(menu.selection) do s
        if s == "ΔR₂"
            connect!(state["heatmap_data"], state["DeltaR2_heatmap"])
            state["heatmap_label"][] = "ΔR₂ (s⁻¹)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["DeltaR2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Relative Intensity"
            connect!(state["heatmap_data"], state["II0_heatmap"])
            state["heatmap_label"][] = "I/I₀"
            state["heatmap_limits"][] = (0., 1.)
            state["heatmap_cm"][] = Reverse(:viridis)
        elseif s == "Reference R₂"
            connect!(state["heatmap_data"], state["refR2_heatmap"])
            state["heatmap_label"][] = "Reference R₂ (s⁻¹)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["refR2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Chemical Shift Perturbation"
            connect!(state["heatmap_data"], state["csps_heatmap"])
            state["heatmap_label"][] = "Chemical Shift Perturbation (ppm)"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["csps_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        elseif s == "Reduced χ²"
            connect!(state["heatmap_data"], state["reducedchi2_heatmap"])
            state["heatmap_label"][] = "Reduced χ²"
            state["heatmap_limits"][] = (0., maximum(filter(!isnan,state["reducedchi2_heatmap"][])))
            state["heatmap_cm"][] = :viridis
        end
    end
    notify(menu.selection)

    ax_heatmap = bottom_panel[2,2] = Axis(fig,
        xlabel="Peak",
        ylabel="Cocktail",
        yreversed=true,
        xzoomlock=true,
        yzoomlock=true,
        xrectzoom=false,
        yrectzoom=false,
        xpanlock=true,
        ypanlock=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax_heatmap, state["heatmap_data"],
        colormap=state["heatmap_cm"],
        colorrange=state["heatmap_limits"])
    cb = Colorbar(bottom_panel[2,3], hm, label=state["heatmap_label"])
    scatter!(ax_heatmap, state["heatmap_point"], color=c4, marker=:star5, markersize=12)
    # Event handler for heatmap click
    on(events(ax_heatmap).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.release
            plt, i = pick(fig)
            plt == hm || return
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
        s = menu.selection[]
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

    colsize!(bottom_panel, 1, Relative(1/4))
    # rowsize!(fig.layout, 1, Relative(1/4))
    rowsize!(fig.layout, 2, Relative(1/3))

    fig
end

