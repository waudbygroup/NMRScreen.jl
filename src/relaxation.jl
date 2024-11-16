function fitrelaxation(spec)
    t = data(spec, TrelaxDim)
    y = data(spec)

    # normalise to first point so I0 ≈ 1
    yerr = ones(length(y)) * spec[:noise] / y[1]
    y = y / y[1]

    # set up model and initial parameters
    model(t, p) = @. p[1] * exp(-p[2] * t)
    p0 = [1., 5.] # I0, R2

    # fitting
    fit = curve_fit(model, t, y, p0)
    I0, R2 = coef(fit)
    I0err, R2err = stderror(fit)

    # store fitted curve for plotting
    tpred = range(0, maximum(t), 100)
    ypred = model(tpred, coef(fit))
    
    # normalise by initial intensity
    ypred = ypred / I0
    y = y / I0
    yerr = yerr / I0
    
    RelaxationResult(R2, R2err, I0, I0err, t, y, yerr, tpred, ypred)
end

Base.show(io::IO, result::RelaxationResult) = print(io, "R2: ", result.R2 ± result.R2_error, " s-1")