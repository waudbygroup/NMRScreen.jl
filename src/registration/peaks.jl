function detectpeaks(spec, snr_threshold=8)
    x = data(spec, F1Dim)
    y = data(spec[:,1]) / spec[:noise]

    pks = findmaxima(y)     # detect peaks
    peakheights!(pks; min=snr_threshold) # filter list

    peak_shifts = x[pks.indices]
    peak_heights = y[pks.indices]
    
    # [peak_shifts peak_heights] # return matrix of peaks
    peak_shifts
end
