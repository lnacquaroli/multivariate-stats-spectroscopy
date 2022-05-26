"""Estimate the percent significance level using Faber's method."""
function fsl(sev, ev, sm, lg)
    Fsl = zeros(sm)
    if lg <= 2000 && sm <= 100
        # lgdf = [ (lg - j + 1) for j=1:sm ]
        lgdf = lg .- collect(1:sm) .+ 1
        # smdf = [ (sm - j + 1) for j=1:sm ]
        smdf = sm .- collect(1:sm) .+ 1
        # fabertbl = fabertable()
        # df1 = [interp2(fabertbl, lgdf[j], smdf[j]) for j=1:sm-1 ]
        df1 = fabertable(lgdf, smdf, sm-1)
        # df2 = [ (lgdf[j]*smdf[j] - df1[j]) for j=1:sm-1 ]
        df2 = lgdf[1:end-1] .* smdf[1:end-1] .- df1
        for j = 1:sm-1
            # df1[j] = interp2(fabertbl,lgdf(k),smdf(k))
            # df2[j] = lgdf[j]*smdf[j] - df1[j]
        # end
        # for j=1:sm-1
            sev = sum(ev[j+1:sm])
            # f = (ev[j]/sev)*(df2[j]/df1[j])
            Fsl[j] = significance_level((ev[j]/sev)*(df2[j]/df1[j]), df1[j], df2[j])
            # if Fsl[j] < 1e-2; Fsl[j] = 0 end
        end
    end
    Fsl[Fsl .< 1.0e-2] .= 0.0
    Fsl[sm] = NaN
    return Fsl
end

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

    # println(lgdf)
    # println(smdf)

    knots = (sort(collect(lgdf[end]:-(lgdf[end]-lgdf[1])/(size(A, 1)-1):lgdf[1])), sort(collect(smdf[end]:-(smdf[end]-smdf[1])/(size(A, 2)-1):smdf[1])))
    itp = interpolate(knots, A, Gridded(Linear()))
    return [ itp(lgdf[j], smdf[j]) for j=1:len_sm ]

end
