% reflectancia y funcion dielectrica IR en un nanocompuesto
% JOURNAL OF RAMAN SPECTROSCOPY
% J. Raman Spectrosc. 2007; 38: 634–646
% Published online 11 April 2007 in Wiley InterScience
% (www.interscience.wiley.com) DOI: 10.1002/jrs.1703

close all
clear all

w = 340:0.1:420;

% modos TO y LO
wto = 367;
wlo = 402;
% damping
g = 1;
% porosidad
f = 0.9;
% constante dielectrica a frecuencia infinita
einf = 8.5; 
%einf = 12; % Si

% modelo de funcion dielectrica
ew = einf*((wlo^2-w.^2-1i*w*g)./(wto^2-w.^2-1i*w*g));

% LLL
%  eeff = looyenga(air(1./w),sqrt(ew),f).^2;

% maxwell-garnett
eeff = ew.*(((2-f)*air(1./w).^2+f*ew)./(f*air(1./w)+(2-f).*ew));

% bruggeman
%eeff = bruggeman(air(1./w),sqrt(ew),f).^2;

figure
plot(w,real(eeff),'r')
hold on
plot(w,imag(eeff),'k')
plot(w,imag(-1./eeff),'b')
ylabel('Funcion dielectrica')
xlabel('Frecuencia')
legend('\Re(\epsilon)','\Im(\epsilon)','\Im(-1/\epsilon)')

% reflectancia a incidencia normal
N = eeff.^0.5;
n = real(N); k = imag(N);
R = ((n-1).^2+k.^2)./((n+1).^2+k.^2);
figure
plot(w,R)
ylabel('Reflectancia')
xlabel('Frecuencia')