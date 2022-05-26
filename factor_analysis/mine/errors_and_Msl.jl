function errors_and_Msl(sev, sdf, ev, sm, lg, rowsX, colsX)
    m = collect(1:sm)
    re = Vector{Float64}(undef, sm)
    ie = similar(re); xe = similar(re); ind = similar(re);
    Msl = similar(re);
    for j = 1:sm-1
        re[j] = sqrt(sev[j+1] / (lg * (sm - j)))
        ie[j] = re[j] * sqrt(j / sm)
        xe[j] = re[j] * sqrt((sm-j) / sm)
        ind[j] = re[j] / (sm-j)^2
        # Estimate the percent significance level using Malinowski's method
        temp1 = (sdf[j+1] * ev[j]) / ((rowsX - j + 1) * (colsX - j + 1) * sev[j+1])
        Msl[j] = significance_level(temp1, 1, sm-j)
    end
    re[end] = NaN; ind[end] = NaN
    Msl[Msl .< 1.0e-2] .= 0; Msl[end] = NaN
    return m, re, ie, xe, ind, Msl
end
