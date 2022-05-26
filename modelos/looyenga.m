function neff = looyenga(n1,n2,p)

% this function calculate the complex refractive index of the porous
% material using Looyenga-Landau-Lifshitz model, for nonmagnetic
% and isotropic materials. Source: IEEE TRANSACTIONS ON GEOSCIENCE
% AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000.
%
% Usage:
%       neff = looyenga(n1,n2,p)
%
% Input:
%       n1:    refractive index of material 1
%       n2:    refractive index of material 2
%       p:     porosity parameter (real number, 0<p<1), proportion of n1
%              over total
%
% Output:
%       neff:  effective refractive index
%
%
% See also maxwellgarnett, monecke, bruggeman

% dielectric function of each media
df1 = n1.^2;
df2 = n2.^2;

% effective medium dielectric, looyenga: IEEE TRANSACTIONS ON GEOSCIENCE
% AND REMOTE SENSING, VOL. 38, NO. 3, MAY 2000
dfeff = ( ( (1-p) * (df2.^(1/3)) ) + ( (df1.^(1/3)) .* p) ).^3;

% compute effective reractive index
neff = dfeff.^0.5;

% EOF looyenga.m
end