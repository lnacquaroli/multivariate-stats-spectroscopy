function [Results] = factor_analysis_pca(d,n)
%factor_analysis_pca.m
% PRINCIPAL FACTOR ANALYSIS - a program designed to help determine
% the number of significant factors in a data matrix.
%
% factor_analysis_pca(d) or factor_analysis_pca(d,n)
%
% d = data matrix
% n = number of principal factors
%
% This program requires uses slf.m and fabertbl.m already (included in the function)
% Main author: Edmund R. Malinowski
%
% Output:
%   Results =  structure with fields:
%       ev = Explained variance
%       re = Real error
%       ie = inbebid_(?) error 
%       xe = xerror
%       ind = Indicator function
%       Msl = Malinowski significance level
%       Fsl = Faber significance level
%       u = Normalised abstract row factors matrix
%       v = Normalised abstract column factors matrix
%       s = Diagonal matrix singular values
%       dr = Reproduced data matrix with principal factors
%       de = Errors data reproduction
%
% modified by: leniac
%

format short e
if nargin == 1, n = 0; end
[r,c] = size(d);
sm = c; lg = r;
if r < c, sm = r; lg = c; end
[u,s,v] = svd(d,0);
for j = 1:sm
    ev(j) = s(j,j) * s(j,j);
    rev(j) = ev(j) / ((r-j+1)*(c-j+1));
end
sev(sm+1) = 0;
sdf(sm+1) = 0;
for k = sm:-1:2
    sev(k) = sev(k+1) + ev(k);
    sdf(k) = sdf(k+1) + (r-k+1) * (c-k+1);
end
for l = 1:sm-1
    m(l) = l;
    re(l) = sqrt(sev(l+1) / (lg * (sm-l)));
    ie(l) = re(l) * sqrt(l/sm);
    xe(l) = re(l) * sqrt((sm-l)/sm);
    ind(l) = re(l) / (sm-l)^2;
end
m(sm) = sm; re(sm) = NaN; ind(sm) = NaN;
out(:,1) = m';
out(:,2) = ev';
out(:,3) = re';
out(:,4) = ind';
% Estimate the percent significance level using Malinowski's method
for j = 1:sm-1
    f = (sdf(j+1) * ev(j)) / ((r-j+1) * (c-j+1) * sev(j+1));
    Msl(j) = slf(f,1,sm-j);
    if Msl(j) < 1e-2, Msl(j) = 0; end
end
Msl(sm) = NaN;
out(:,5) = Msl';
% Estimate the percent significance level using Faber's method
for k = 1:sm; Fsl(k) = NaN; end
if lg <= 2000 && sm <= 100
    for i=1:sm
        lgdf(i) = lg-i+1;
        smdf(i) = sm-i+1;
    end
    fabertbl = fabertable(); %load fabertbl
    for k=1:sm-1
        %df1(k) = table2(fabertbl,lgdf(k),smdf(k));
        df1(k) = interp2(fabertbl,lgdf(k),smdf(k));
        df2(k) = lgdf(k)*smdf(k) - df1(k);
    end
    for j=1:sm-1
        sev = sum(ev(j+1:sm));
        f(j) = (ev(j)/sev)*(df2(j)/df1(j));
        Fsl(j) = slf(f(j),df1(j),df2(j));
        if Fsl(j) < 1e-2, Fsl(j) = 0; end
    end
end
Fsl(sm) = NaN;
% Output
out(:,6) = Fsl';
clc
disp('                          Principal Factor Analysis')
disp(' ')
disp('       n           EV           RE           IND    %SL(Malinowski) %SL(Faber)')
disp(out)
if n > 0
    % Results
    u = u(:,1:n);
    s = s(1:n,1:n);
    v = v(:,1:n);
    dr = u * s * v';
    de = dr - d;
    % Save errors and results to return
    Results = struct('ev', ev, 're', re, 'ie', ie, 'xe', xe, 'ind', ind,...
                     'Msl', Msl, 'Fsl', Fsl, 'u', u, 'v', v, 's', s, 'dr', dr, 'de', de);
else
    Results = [];
end
end

function [sl] = slf(f,df1,df2)
% slf.m - SIGNIFICANCE LEVELS FOR F-RATIOS - a program designed to calculate
%         the percent significance level for a specified F(df1,df2).
% @ copyright 2000 by Edmund R. Malinowski
% f = fisher F value
% df1 = numerator degrees of freedom
% df2 = denominator degrees of freedom
format short e
if f > 1e5, f = 1e5; end
df(1) = df1;
df(2) = df2;
df(3) = df(1) + df(2);
xx = (df(1) * f) / (df(2) + df(1) * f);
d1 = 0.5 * df(1);
d2 = 0.5 * df(2);
d3 = d1 + d2;
cl = xx;
cx = 1 - xx;
if d1 >= (d3 * xx)
    xxx = xx;
    yy = d1;
    zz = d2;
    invert = 1;
else
    xxx = cx;
    cx = xx;
    yy = d2;
    zz = d1;
    invert = 0;
end
tt = 1;
aa = 1;
cl = 1;
nn = floor(zz + cx * d3);
rr = xxx / cx;
while nn >= 0
    tmp = zz - aa;
    if nn == 0, rr = xxx; end
    burp = 1;
    while burp == 1
        tt = tt * tmp * rr / (yy + aa);
        cl = cl + tt;
        tmp = abs(tt);
        skip = 0;
        if (tmp <= 1e-8), skip = 1; end
        if (tmp <= cl * 1e-8), skip = skip + 1; end
        if skip == 2, burp = 0; nn = 0; end
        aa = aa + 1;
        nn = nn - 1;
        if nn < 0, tmp = d3; d3 = d3 + 1; end
        if nn >= 0, burp = 0; end
    end
end
g = zeros(3);
for i = 1:3
    dxx = 0.5 * df(i);
    if dxx > 0
        af = 0;
        if dxx < 7
            af = 1;
            az = dxx - 1;
            az = az + 1;
            while az < 7
                dxx = az;
                af = af * az;
                az = az + 1;
            end
            dxx = dxx +1;
            af = -log(af);
        end
        az = 1 / dxx^2;
        g(i) = af + (dxx - 0.5) * log(dxx) - dxx + 0.918939;
        g(i) = g(i) + (((-0.000595238*az+0.00079365)*az-0.00277778)*az+0.0833333)/dxx;
    end
    if dxx <= 0, i = 4; end
end
ggg = g(1) + g(2) - g(3);
cl = cl * exp(yy * log(xxx) + (zz - 1) * log(cx) - ggg) / yy;
if invert == 0, cl = 1 - cl; end
sl = 100 * (1 - cl);
end


function [fbtbl] = fabertable()
fbtbl = [0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,80,100;
    1,1,2,3,4,5,6,7,8,9,10,20,30,40,50,80,100;
    2,2,4,5,7,8,9,10,11,13,14,26,36,47,59,92,113;
    3,3,5,7,8,10,11,13,14,15,17,29,42,54,66,98,121;
    4,4,7,8,10,12,13,15,16,18,19,33,46,59,71,105,128;
    5,5,8,10,12,14,15,17,19,20,22,36,49,63,75,111,134;
    6,6,9,11,13,15,17,19,21,22,24,39,53,66,79,116,140;
    7,7,10,13,15,17,19,21,23,25,26,42,56,69,83,121,145;
    8,8,11,14,16,19,21,23,25,27,28,44,59,73,87,125,150;
    9,9,13,15,18,20,22,25,27,29,30,47,62,76,90,129,154;
    10,10,14,17,19,22,24,26,28,30,32,50,65,79,93,133,159;
    20,20,26,29,33,36,39,42,44,47,50,71,88,105,121,167,195;
    30,30,36,42,46,49,53,56,59,62,65,88,109,128,146,195,226;
    40,40,47,54,59,62,66,69,73,76,79,105,128,147,168,219,252;
    50,50,59,66,71,75,79,83,87,90,93,121,146,168,187,242,276;
    60,60,70,77,82,87,92,95,100,104,107,137,162,186,206,265,300;
    70,70,80,88,93,100,104,109,112,116,120,152,179,202,226,286,322;
    80,80,92,98,105,111,116,121,125,129,133,167,195,219,242,304,344;
    90,90,101,110,117,123,128,132,138,142,146,182,211,236,260,324,363;
    100,100,113,121,128,134,140,145,150,155,159,195,226,252,276,344,384;
    200,200,218,227,240,248,255,262,269,275,281,329,368,402,432,515,565;
    300,300,320,338,348,358,368,377,383,391,398,455,500,541,576,670,726;
    400,400,425,442,456,467,478,488,495,504,513,576,628,670,712,817,878;
    500,500,527,546,562,573,586,596,607,616,625,695,752,799,844,957,1024;
    600,600,630,651,667,681,693,707,716,728,737,813,874,924,972,1095,1167
    700,700,734,756,774,789,801,814,825,837,846,927,991,1048,1098,1229,1305;
    800,800,836,860,878,893,908,922,934,947,954,1044,1111,1170,1222,1359,1440;
    900,900,937,964,984,1002,1017,1028,1043,1055,1064,1155,1229,1291,1345,1490,1573;
    1000,1000,1038,1066,1087,1106,1119,1137,1149,1161,1174,1270,1345,1410,1468,1617,1704;
    1500,1500,1550,1582,1607,1629,1646,1667,1680,1697,1710,1828,1917,1993,2064,2240,2342;
    2000,2000,2055,2091,2123,2150,2150,2170,2191,2227,2242,2376,2478,2566,2644,2842,2958];

end
