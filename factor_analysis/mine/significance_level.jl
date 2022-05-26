function significance_level(f, df1, df2)
    # SIGNIFICANCE LEVELS FOR F-RATIOS - a program designed to calculate the percent significance level for a specified F(df1,df2).
    # @ copyright 2000 by Edmund R. Malinowski
    # f = fisher F value
    # df1 = numerator degrees of freedom
    # df2 = denominator degrees of freedom
    if f > 1e5; f = 1e5 end
    df = [df1; df2; df1+df2]
    xx = (df[1] * f) / (df[2] + df[1] * f)
    d1 = 0.5 * df[1]; d2 = 0.5 * df[2]; d3 = d1 + d2
    cl = xx; cx = 1 - xx
    if d1 >= (d3 * xx)
        xxx = xx
        yy = d1
        zz = d2
        invert = 1
    else
        xxx = cx
        cx = xx
        yy = d2
        zz = d1
        invert = 0
    end
    tt = 1; aa = 1; cl = 1
    nn = floor(zz + cx * d3)
    rr = xxx / cx
    while nn >= 0
        tmp = zz - aa
        if nn == 0; rr = xxx end
        burp = 1
        while burp == 1
            tt = tt * tmp * rr / (yy + aa)
            cl = cl + tt
            tmp = abs(tt)
            skip = 0
            if tmp <= 1e-8; skip = 1 end
            if (tmp <= cl * 1e-8); skip = skip + 1 end
            if skip == 2; burp = 0; nn = 0 end
            aa = aa + 1
            nn = nn - 1
            if nn < 0; tmp = d3; d3 = d3 + 1 end
            if nn >= 0; burp = 0 end
        end
    end
    g = zeros(3)
    for i = 1:3
        dxx = 0.5 * df[i]
        if dxx > 0
            af = 0
            if dxx < 7
                af = 1
                az = dxx - 1
                az = az + 1
                while az < 7
                    dxx = az
                    af = af * az
                    az = az + 1
                end
                dxx = dxx + 1
                af = -log(af)
            end
            az = 1 / dxx^2
            g[i] = af + (dxx - 0.5) * log(dxx) - dxx + 0.918939
            g[i] = g[i] + (((-0.000595238*az + 0.00079365)*az - 0.00277778)*az + 0.0833333)/dxx
        end
        if dxx <= 0; i = 4 end
    end
    ggg = g[1] + g[2] - g[3]
    cl = cl * exp(yy * log(xxx) + (zz - 1) * log(cx) - ggg) / yy
    if invert == 0; cl = 1 - cl end

    return 100*(1 - cl)

end #significance_level
