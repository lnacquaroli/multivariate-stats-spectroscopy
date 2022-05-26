close all
clear
clc

% ingresar la matriz de errores
Merror = load('MatrizError.dat');

% ingresar la matriz de datos
Mdatos = load('MatrizDatos.dat');

% datos del tiempo
tiempo = load('ti2.dat');

% x-axis
nu = 1:size(Merror,1);

% aplica factor principal a la matriz de errores
pfa2(Merror);

% ingresar la matriz de datos
% Mdatos = load('pru3.dat');
% tiempo =1:length(Mdatos(1,:));
% aplica factor principal a la matris de datos
pfa2(Mdatos);

% ingresar numero de factores
n = 2;

% aplica factor principal a la matris de datos con n factores
pfa2(Mdatos, n);

% carga resultados temporales
load temp
load errors
figure();
plot(dr)
title('Reproduccion de datos abstractos')

%anfa4a;

% Numero de espectro para los que hay targets
nes = [1; size(Mdatos,2)];

% Concentation inicial de los espectros seleccionados
ci = [1; 0];
cf = 1 - sum(ci, 2);
b = [ci cf];
b = pinv(b);
spectra = Mdatos(:, nes);
xx = spectra * b;

% calcula la matriz fila
R = u * s;

% calcula la matriz diagonal
Mdiagonal = diag(ev(1:n).^(-1));
T = Mdiagonal * R' * xx;

% calcula la matriz de factores reales
X = R * T;

% calcula la matriz de concentraciones
Y = T^(-1) * v';


tiempo = tiempo - tiempo(1);

figure()
plot(nu,X,'-')
ylabel('Componentes principales reales')
xlabel('Numero de onda [cm^{-1}]')

figure();
plot(tiempo, Y', 'o')
ylabel('Concentraciones principales reales')
xlabel('Tiempo')
RRR = X * Y;

figure()
plot(nu,RRR,'-')
% ylabel('Componentes principales reales')
xlabel('Tiempo')
title('Reproduccion de datos reales')
% RMS = norm(dr-RRR)/norm(dr)
