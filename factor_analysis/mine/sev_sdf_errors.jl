function sev_sdf_errors(ev, sm, rowsX, colsX)
    sev = zeros(sm + 1); sdf = zeros(sm + 1)
    for j = sm:-1:2
        sev[j] = sev[j+1] + ev[j]
        sdf[j] = sdf[j+1] + (rowsX - j + 1)*(colsX - j + 1)
    end
    return sev, sdf
end
