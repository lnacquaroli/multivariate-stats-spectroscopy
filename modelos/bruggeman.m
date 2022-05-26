function neff = bruggeman(n1,n2,p)

% this function calculate the complex refractive index of a porous material
% using Bruggeman model, for nonmagnetic and isotropic materials. Source:
% PHYSICAL REVIEW B VOLUME 61, NUMBER 15 15 APRIL 2000-I.
%
% Usage:
%       neff = bruggeman(n1,n2,p)
%
% Input:
%       n1:   refractive index of embedded material
%       n2:   refractive index of host material
%       p:    porosity parameter (real number, 0<p<1)
%
% Output:
%       neff:  effective refractive index
%
% See also looyenga, monecke, maxwellgarnett

% dielectric function of each media
df1 = n1.^2;
df2 = n2.^2;

% effective dielectric function: PHYSICAL REVIEW B VOLUME 61, NUMBER 15 15 APRIL 2000-I. solved with mathematica for dfeff.
aux1 = df1 ./ df2 ;
aux2 = 0.5 .* ( ( ( (1.5 .* p) - 0.5 ) .* (1-aux1) ) + ( 0.5 .* aux1 ) ) ;
aux3 = aux2 + ( aux2.^2 + (0.5 .* aux1) ).^0.5 ;
dfeff = df2 .* aux3;

% compute effective refractive index
neff = dfeff.^0.5;

% EOF
end