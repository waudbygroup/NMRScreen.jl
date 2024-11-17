function makescreeningpeaks(refspec, boundspec, cocktail_id, librarypeak_ids, librarypeaks, refpeaks, boundpeaks)
    # 1. match reference peaks to library
    x_offset, y_scale, matches, unmatched_expected, unmatched_observed = rpm(librarypeaks, refpeaks; intensity_weight=0.)
    # matches = pairs of (library, ref) peaks
    # unmatched_expected = library peaks without a match
    # unmatched_observed = ref peaks without a match
    matched_refpeaks = refpeaks[[m[2] for m in matches], :]
    matched_librarypeaks = librarypeaks[[m[1] for m in matches], :]
    matched_ids = librarypeak_ids[[m[1] for m in matches]]
    # TODO include unrecognised peaks from refpeaks and give them 'unknown' labels
    missing_librarypeak_ids = librarypeak_ids[unmatched_expected]
    @warn missing_librarypeak_ids

    # 2. match bound peaks to reference
    x_offset2, y_scale2, matches2, unmatched_expected2, unmatched_observed2 = rpm(matched_refpeaks, boundpeaks; intensity_weight=0.000)
    matched_boundpeaks = boundpeaks[[m[2] for m in matches2], 1] |> vec
    matched_refpeaks2 = matched_refpeaks[[m[1] for m in matches2], 1] |> vec
    matched_librarypeaks2 = matched_librarypeaks[[m[1] for m in matches2], 1] |> vec
    matched_ids2 = matched_ids[[m[1] for m in matches2]]
    missing_refpeak_ids = matched_ids[unmatched_expected2]
    @warn missing_refpeak_ids

    # 3. create ScreeningPeak objects
    screeningpeaks = []
    for (i, peak_id) in enumerate(matched_ids2)
        libpeak = matched_librarypeaks2[i]
        refpeak = matched_refpeaks2[i]
        boundpeak = matched_boundpeaks[i]
        newpeak = ScreeningPeak(refspec, boundspec, peak_id, peak_id[1:end-1], cocktail_id, libpeak, refpeak, boundpeak)
        push!(screeningpeaks, newpeak)
    end
    reverse!(screeningpeaks)

    return screeningpeaks, missing_librarypeak_ids
end


function ScreeningPeak(refspec,
        boundspec,
        peak_id,
        fragment_id,
        cocktail_id,
        library_shift,
        reference_shift,
        bound_shift)

    reference_intensity = refspec[Near(reference_shift), 1] / scale(refspec)
    reference_intensity_error = refspec[:noise] / scale(refspec)
    bound_intensity = boundspec[Near(bound_shift), 1] / scale(boundspec)
    bound_intensity_error = boundspec[:noise] / scale(boundspec)

    reference_relaxation = fitrelaxation(refspec[Near(reference_shift),:])
    bound_relaxation = fitrelaxation(boundspec[Near(bound_shift),:])

    ScreeningPeak(refspec,
        boundspec,
        peak_id,
        fragment_id,
        cocktail_id,
        library_shift,
        reference_shift,
        bound_shift,
        reference_intensity,
        reference_intensity_error,
        bound_intensity,
        bound_intensity_error,
        reference_relaxation,
        bound_relaxation)
end

Base.show(io::IO, sp::ScreeningPeak) = begin
    print(io, "ScreeningPeak(", sp.cocktail_id, "/", sp.peak_id, 
          ", ", round(sp.reference_shift,digits=2), 
          "ppm, II0: ", II0(sp), ", ΔR2: ", DeltaR2(sp), " s-1)")
end

function II0(sp::ScreeningPeak)
    I = sp.bound_intensity ± sp.bound_intensity_error
    I0 = sp.reference_intensity ± sp.reference_intensity_error
    I / I0
end

function refR2(sp::ScreeningPeak)
    sp.reference_relaxation.R2 ± sp.reference_relaxation.R2_error
end

function DeltaR2(sp::ScreeningPeak)
    R2bound = sp.bound_relaxation.R2 ± sp.bound_relaxation.R2_error
    R2free = sp.reference_relaxation.R2 ± sp.reference_relaxation.R2_error
    R2bound - R2free
end

function csp(sp::ScreeningPeak)
    sp.bound_shift - sp.reference_shift
end

function reducedchi2(sp::ScreeningPeak)
    (sp.reference_relaxation.reducedchi2 + sp.bound_relaxation.reducedchi2) / 2
end