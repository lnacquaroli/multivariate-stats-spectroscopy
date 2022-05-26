module PFAMalinowski

export PFA, PFATT

using LinearAlgebra
using Interpolations
using PrettyTables

"""Wrap the output into different types."""
abstract type Output end
struct ErrorsOutput{T0,T1} <: Output where {T0<:Int64, T1<:Float64}
    index::Array{T0}; explvar::Array{T1}; realerr::Array{T1}; indfunc::Array{T1};
    Msl::Array{T1}; Fsl::Array{T1};
end
struct MatrixOutput{T1} <: Output where {T1<:Float64}
    U::Array{T1}; S::Array{T1}; V::Array{T1}; Vt::Array{T1};
    Xreprod::Array{T1}; Xerror::Array{T1};
end
struct PFAResults{T1, T2} <:Output where {T1<:ErrorsOutput, T2<:MatrixOutput}
    errorsResults::T1; matrixResults::T2;
end
struct PFATTResults{T1} <:Output where {T1<:Float64}
    row::Array{T1}; diag::Array{T1}; transf::Array{T1}; realFac::Array{T1};
    realConc::Array{T1}; realSpectraReprod::Array{T1}; rms::T1;
end

"""PRINCIPAL FACTOR ANALYSIS
Determine the number of significant factors in a data matrix.
X: cols = observations (samples), rows = variables (wavelengths, etc)."""
function PFA(X::AbstractArray, k::Int64=0)
    rowsX, colsX = size(X)
    sm = colsX; lg = rowsX
    if rowsX < colsX; sm = rowsX; lg = colsX end
    # Perform thin-SVD
    F, ev, rev = ev_rev_svd(X, rowsX, colsX, sm)
    # Errors and Malinowski significance level
    m, re, ie, xe, ind, Msl = errors_Msl(ev, sm, lg, rowsX, colsX)
    # Faber significance level
    Fsl = fsl(ev, sm, lg)
    # Print results
    results_errors_out = ErrorsOutput(m, ev, re, ind, Msl, Fsl)
    print_results(results_errors_out)
    # Return results
    return k > 0 ? PFAResults(results_errors_out, MatrixOutput(F.U[:,1:k], Matrix(Diagonal(F.S))[1:k,1:k], F.V[:,1:k], F.Vt[:,1:k], F.U[:,1:k] * Matrix(Diagonal(F.S))[1:k,1:k] * F.V[:,1:k]', (F.U[:,1:k] * Matrix(Diagonal(F.S))[1:k,1:k] * F.V[:,1:k]') .- X)) : []
end #PFA

"""Perform target transformation."""
function PFATT(spectra_target, initial_conc, X, A, k)
    # Selected spectra
    # spectra_target = spectra_target[:]
    spectra = X[:, spectra_target[:]]
    # Estimate final concentration
    initial_conc = initial_conc[:]
    final_conc = 1.0 .- sum(initial_conc, dims=2)
    # Concentration matrix
    conc_matrix = pinv([initial_conc final_conc])
    # conc_matrix = pinv(conc_matrix)
    # Row matrix
    Row_matrix = A.matrixResults.U * A.matrixResults.S
    # Diagonal matrix
    Diag_matrix = Matrix(Diagonal(inv.(A.errorsResults.explvar[1:k])))
    # Transformation matrix
    tmp1 = spectra * conc_matrix
    Transformation_matrix = Diag_matrix * Row_matrix' * tmp1
    # Real factor matrix
    Real_factors_matrix = Row_matrix * Transformation_matrix
    # Real concentration matrix
    Real_conc_matrix = inv(Transformation_matrix) * A.matrixResults.V'
    # Real spectra reproduction
    Real_spectra_reprod = Real_factors_matrix * Real_conc_matrix
    # Error in reproduction
    rms = norm(A.matrixResults.Xreprod .- Real_spectra_reprod) / norm(A.matrixResults.Xreprod)
    # Wrap up results
    PFATTResults(Row_matrix, Diag_matrix, Transformation_matrix, Real_factors_matrix,
                 Real_conc_matrix, Real_spectra_reprod, rms)
end #PFATT

"""Display table results from PFA."""
function print_results(errors)
    header = ["Factors" "EV" "RE" "IND" "%SL (Malinowski)" "%SL (Faber)"]
    hc = crayon"bold blue"
    n = errors.index
    ev = errors.explvar
    re = errors.realerr
    ind = errors.indfunc
    msl = round.(errors.Msl, digits=4)
    fsl = round.(errors.Fsl, digits=4)
    data = [n ev re ind msl fsl]
    h1 = Highlighter(f = (data,i,j) -> (isnan(data[i,j])), crayon=crayon"bold red")
    pretty_table(data, header;
                 alignment=:c, header_crayon=hc, formatter=ft_printf("%0.4e", [2,3,4]),
                 highlighters=(h1), crop=:none)
end #print_results

"""Perform thin-SVD."""
function ev_rev_svd(X, rowsX, colsX, sm)
    F = svd(X, full=false) # F.U, F.S, F.V
    ev = (F.S).^2
    temp1 = collect(1:sm)
    rev = ev ./ ((rowsX .- temp1 .+ 1 ) .* (colsX .- temp1 .+ 1))
    return F, ev, rev
end #ev_rev_svd

"""Errors and Malinowski significance level."""
function errors_Msl(ev, sm, lg, rowsX, colsX)
    sev = zeros(sm + 1); sdf = zeros(sm + 1)
    for j = sm:-1:2
        sev[j] = sev[j+1] + ev[j]
        sdf[j] = sdf[j+1] + (rowsX - j + 1)*(colsX - j + 1)
    end
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
end #errors_Msl

"""SIGNIFICANCE LEVELS FOR F-RATIOS.
Careful though: Scientists rise up against statistical significance https://www.nature.com/articles/d41586-019-00857-9
"""
function significance_level(f, df1, df2)
    # Calculate the percent significance level for a specified F(df1,df2).
    # by Edmund R. Malinowski
    # f = fisher F value
    # df1 = numerator degrees of freedom
    # df2 = denominator degrees of freedom
    if f > 1e5; f = 1e5 end
    df = [df1; df2; df1+df2]
    xx = (df[1] * f) / (df[2] + df[1] * f)
    d1 = 0.5 * df[1]; d2 = 0.5 * df[2]; d3 = d1 + d2
    cl = xx; cx = 1 - xx
    if d1 >= (d3 * xx)
        xxx = xx; yy = d1; zz = d2
        invert = 1
    else
        xxx = cx; cx = xx; yy = d2; zz = d1
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
            aa = aa + 1; nn = nn - 1
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
                af = 1; az = dxx - 1; az = az + 1
                while az < 7
                    dxx = az; af = af * az; az = az + 1
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

"""Estimate the percent significance level using Faber's method."""
function fsl(ev, sm, lg)
    Fsl = zeros(sm)
    if lg <= 2000 && sm <= 100
        temp1 = collect(1:sm)
        lgdf = lg .- temp1 .+ 1
        smdf = sm .- temp1 .+ 1
        df1 = fabertable(lgdf, smdf, sm-1)
        df2 = lgdf[1:end-1] .* smdf[1:end-1] .- df1
        for j = 1:sm-1
            sev = sum(ev[j+1:sm])
            Fsl[j] = significance_level((ev[j]/sev)*(df2[j]/df1[j]), df1[j], df2[j])
        end
    end
    Fsl[Fsl .< 1.0e-2] .= 0.0
    Fsl[sm] = NaN
    return Fsl
end #fsl

"""Interpolate the degree of freedom obtained from several Monte Carlo simulations.
N. Faber, J. Chemometrics 2000, 14: 371-374 and E. Malinowski, J. Chemometrics 13, 69-81 (1999)."""
function fabertable(lgdf, smdf, len_sm)
    A = [0 1 2 3 4 5 6 7 8 9 10 20 30 40 50 80 100;
    1 1 2 3 4 5 6 7 8 9 10 20 30 40 50 80 100;
    2 2 4 5 7 8 9 10 11 13 14 26 36 47 59 92 113;
    3 3 5 7 8 10 11 13 14 15 17 29 42 54 66 98 121;
    4 4 7 8 10 12 13 15 16 18 19 33 46 59 71 105 128;
    5 5 8 10 12 14 15 17 19 20 22 36 49 63 75 111 134;
    6 6 9 11 13 15 17 19 21 22 24 39 53 66 79 116 140;
    7 7 10 13 15 17 19 21 23 25 26 42 56 69 83 121 145;
    8 8 11 14 16 19 21 23 25 27 28 44 59 73 87 125 150;
    9 9 13 15 18 20 22 25 27 29 30 47 62 76 90 129 154;
    10 10 14 17 19 22 24 26 28 30 32 50 65 79 93 133 159;
    20 20 26 29 33 36 39 42 44 47 50 71 88 105 121 167 195;
    30 30 36 42 46 49 53 56 59 62 65 88 109 128 146 195 226;
    40 40 47 54 59 62 66 69 73 76 79 105 128 147 168 219 252;
    50 50 59 66 71 75 79 83 87 90 93 121 146 168 187 242 276;
    60 60 70 77 82 87 92 95 100 104 107 137 162 186 206 265 300;
    70 70 80 88 93 100 104 109 112 116 120 152 179 202 226 286 322;
    80 80 92 98 105 111 116 121 125 129 133 167 195 219 242 304 344;
    90 90 101 110 117 123 128 132 138 142 146 182 211 236 260 324 363;
    100 100 113 121 128 134 140 145 150 155 159 195 226 252 276 344 384;
    200 200 218 227 240 248 255 262 269 275 281 329 368 402 432 515 565;
    300 300 320 338 348 358 368 377 383 391 398 455 500 541 576 670 726;
    400 400 425 442 456 467 478 488 495 504 513 576 628 670 712 817 878;
    500 500 527 546 562 573 586 596 607 616 625 695 752 799 844 957 1024;
    600 600 630 651 667 681 693 707 716 728 737 813 874 924 972 1095 1167
    700 700 734 756 774 789 801 814 825 837 846 927 991 1048 1098 1229 1305;
    800 800 836 860 878 893 908 922 934 947 954 1044 1111 1170 1222 1359 1440;
    900 900 937 964 984 1002 1017 1028 1043 1055 1064 1155 1229 1291 1345 1490 1573;
    1000 1000 1038 1066 1087 1106 1119 1137 1149 1161 1174 1270 1345 1410 1468 1617 1704;
    1500 1500 1550 1582 1607 1629 1646 1667 1680 1697 1710 1828 1917 1993 2064 2240 2342;
    2000 2000 2055 2091 2123 2150 2150 2170 2191 2227 2242 2376 2478 2566 2644 2842 2958]

    knots = (sort(collect(lgdf[end]:-(lgdf[end]-lgdf[1])/(size(A, 1)-1):lgdf[1])), sort(collect(smdf[end]:-(smdf[end]-smdf[1])/(size(A, 2)-1):smdf[1])))
    itp = interpolate(knots, A, Gridded(Linear()))
    return [ itp(lgdf[j], smdf[j]) for j=1:len_sm ]

end #fabertable

end #PFAMalinowski
