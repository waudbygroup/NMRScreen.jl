function prepcocktails(cocktails)
    # input: cocktails where peaks are BasicPeaks
    # output: cocktails where peaks are ScreeningPeaks
    for cocktail in cocktails
        for (i, peak) in enumerate(cocktail.peaks)
            newpeak = ScreeningPeak(peak)
            cocktail.peaks[i] = newpeak
        end
        sort!(cocktail.peaks, by=sp -> sp.reference_shift, rev=true)
    end

    return cocktails
end