"""
Functions for saving NMR screening data
"""

function write_results(state)
    od = state["config"].output_directory
    if isdir(od)
        backup_od = od * "_previous"
        if isdir(backup_od)
            rm(backup_od; force=true, recursive=true)
        end
        mv(od, backup_od)
    end
    mkdir(od)

    # bgcolor = state["fig"].scene.backgroundcolor[]
    # state["fig"].scene.backgroundcolor[] = RGBf(1,.85,.85)
    write_results_csv(joinpath(od, "results.csv"), state["cocktails"])
    write_heatmaps(state)
    write_peaks(state)
    write_top_fragments(state)

    @info "Finished save"

    # state["fig"].scene.backgroundcolor[] = bgcolor
end

"""
    write_results(filename, cocktails)

Write screening results to a CSV file.
"""
function write_results_csv(filename, cocktails)
    @info "Writing results to $filename"
    open(filename, "w") do io
        # Write the header
        println(io, "cocktail_id,fragment_id,peak_id,library_shift,reference_shift,bound_shift,intensity_ratio,intensity_ratio_error,chemical_shift_perturbation,reference_R2,reference_R2_error,bound_R2,bound_R2_error,DeltaR2,DeltaR2_error,reduced_chi2")
        
        for cocktail in cocktails
            for peak in cocktail.peaks
                println(io, join([
                    peak.cocktail_id,
                    peak.fragment_id,
                    peak.peak_id,
                    peak.library_shift,
                    peak.reference_shift,
                    peak.bound_shift,
                    Measurements.value(II0(peak)),
                    Measurements.uncertainty(II0(peak)),
                    csp(peak),
                    peak.reference_relaxation.R2,
                    peak.reference_relaxation.R2_error,
                    peak.bound_relaxation.R2,
                    peak.bound_relaxation.R2_error,
                    Measurements.value(DeltaR2(peak)),
                    Measurements.uncertainty(DeltaR2(peak)),
                    reducedchi2(peak)
                ], ","))
            end
        end
    end
    return nothing
end


function write_heatmaps(state)
    od = state["config"].output_directory
    od = joinpath(od, "heatmaps")
    if !isdir(od)
        mkdir(od)
    end

    # Switch to CairoMakie
    CairoMakie.activate!()
    

    # DeltaR2
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel="Peak",
        ylabel="Cocktail",
        backgroundcolor=:grey10,
        yreversed=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax, state["DeltaR2_heatmap"],
        colormap=:viridis,
        colorrange=(0., maximum(filter(!isnan,state["DeltaR2_heatmap"][]))))
    cb = Colorbar(fig[1,2], hm, label="ΔR₂ (s⁻¹)")

    # Save the figure
    filename = joinpath(od, "DeltaR2.pdf")
    @info "Saving ΔR2 heatmap to $filename"
    Makie.save(filename, fig)


    # Reference R2
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel="Peak",
        ylabel="Cocktail",
        backgroundcolor=:grey10,
        yreversed=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax, state["refR2_heatmap"],
        colormap=:viridis,
        colorrange=(0., maximum(filter(!isnan,state["refR2_heatmap"][]))))
    cb = Colorbar(fig[1,2], hm, label="ΔR₂ (s⁻¹)")

    # Save the figure
    filename = joinpath(od, "refR2.pdf")
    @info "Saving reference R2 heatmap to $filename"
    Makie.save(filename, fig)


    # II0
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel="Peak",
        ylabel="Cocktail",
        backgroundcolor=:grey10,
        yreversed=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax, state["II0_heatmap"],
        colormap=Reverse(:viridis),
        colorrange=(0., 1))
    cb = Colorbar(fig[1,2], hm, label="Relative intensity")

    # Save the figure
    filename = joinpath(od, "II0.pdf")
    @info "Saving I/I₀ heatmap to $filename"
    Makie.save(filename, fig)


    # CSPs
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel="Peak",
        ylabel="Cocktail",
        backgroundcolor=:grey10,
        yreversed=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax, state["csps_heatmap"],
        colormap=:viridis,
        colorrange=(0., maximum(filter(!isnan,state["csps_heatmap"][]))))
    cb = Colorbar(fig[1,2], hm, label="Chemical shift perturbation (ppm)")

    # Save the figure
    filename = joinpath(od, "csps.pdf")
    @info "Saving csps heatmap to $filename"
    Makie.save(filename, fig)



    # reduced chi2
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel="Peak",
        ylabel="Cocktail",
        backgroundcolor=:grey10,
        yreversed=true,
        yticks=(1:state["n_cocktails"], state["cocktail_ids"])
        )
    hm = heatmap!(ax, state["reducedchi2_heatmap"],
        colormap=:viridis,
        colorrange=(0., maximum(filter(!isnan,state["reducedchi2_heatmap"][]))))
    cb = Colorbar(fig[1,2], hm, label="Reduced χ²")

    # Save the figure
    filename = joinpath(od, "reducedchi2.pdf")
    @info "Saving reduced chi2 heatmap to $filename"
    Makie.save(filename, fig)

    # Switch back to GLMakie
    GLMakie.activate!()
end


function write_peaks(state)
    for cocktail in state["cocktails"]
        for peak in cocktail.peaks
            write_peak(peak, state)
        end
    end
end


function write_peak(peak, state)
    od = state["config"].output_directory
    od = joinpath(od, "peaks")
    if !isdir(od)
        mkdir(od)
    end

    peak_id = peak.peak_id

    # Switch to CairoMakie
    CairoMakie.activate!()
    
    fig = Figure(size=(600,400))
    ax = Axis(fig[1:2,1],
        xlabel="Relaxation time (ms)",
        ylabel="Intensity",
        limits=(nothing, (-.1, 1.1)))
    lines!(ax, peak.reference_relaxation.typred)
    lines!(ax, peak.bound_relaxation.typred)
    errorbars!(ax, peak.reference_relaxation.tye)
    errorbars!(ax, peak.bound_relaxation.tye)
    scatter!(ax, peak.reference_relaxation.ty, label = "Reference")
    scatter!(ax, peak.bound_relaxation.ty, label = "Bound")
    axislegend(ax)

    # info panel
    infotext = """Cocktail: $(peak.cocktail_id)
        Peak: $(peak.peak_id)
        
        ΔR₂: $(DeltaR2(peak)) s⁻¹ 
        I/I₀: $(II0(peak))
        Δδ: $(round(csp(peak), digits=2)) ppm
        Ref. R₂: $(refR2(peak)) s⁻¹
        """
    info = Label(fig[1,2], text = infotext, justification=:left, lineheight=1.2)

    # structure
    smiles = state["library"].fragments[peak.fragment_id].smiles
    img = smilestoimage(smiles)
    ax_structure = Axis(fig[2,2], aspect = DataAspect(), backgroundcolor = :white)
    hidespines!(ax_structure)
    hidedecorations!(ax_structure)
    structure = image!(ax_structure, img)

    # Save the figure
    filename = joinpath(od, "$peak_id.pdf")
    @info "Saving peak output to $filename"
    Makie.save(filename, fig)

    # Switch back to GLMakie
    GLMakie.activate!()
end


function write_top_fragments(state, n=15)
    od = state["config"].output_directory
    
    # Extract top n fragments by DeltaR2
    peaks = [peak for cocktail in state["cocktails"] for peak in cocktail.peaks]
    sorted_peaks = sort(peaks, by=(i->Measurements.value(DeltaR2(i))), rev=true)
    top_peaks = sorted_peaks[1:n]

    # Switch to CairoMakie
    CairoMakie.activate!()

    fig = Figure(size=(800,400))
    ax = Axis(fig[1,1],
        xlabel="ΔR₂ (s⁻¹)",
        ylabel="Peak ID",
        yreversed=true,
        title="Top $n fragments by ΔR₂",
        yticks=(1:n, [peak.peak_id for peak in top_peaks])
    )

    # Plot horizontal bars with error bars
    barpositions = 1:n
    barvalues = [Measurements.value(DeltaR2(peak)) for peak in top_peaks]
    barerrors = [Measurements.uncertainty(DeltaR2(peak)) for peak in top_peaks]
    errorbars!(ax, barvalues, barpositions, barerrors, direction=:x, color=:black, whiskerwidth=5)
    barplot!(ax, barpositions, barvalues, direction=:x)#, color=Makie.wong_colors()[1])

    ax2 = Axis(fig[1,2],
        xlabel="ΔR₂ (s⁻¹)",
        ylabel="Frequency",
        title="Distribution of ΔR₂ values"
    )
    hist!(ax2, [Measurements.value(DeltaR2(peak)) for peak in peaks], bins=20)

    # Save the figure
    filename = joinpath(od, "top_fragments.pdf")
    @info "Saving top fragments bar plot to $filename"
    Makie.save(filename, fig)

    # Switch back to GLMakie
    GLMakie.activate!()
end