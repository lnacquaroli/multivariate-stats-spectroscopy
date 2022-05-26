close all
clear
clc

tmp1 = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/datasetSi_n_with_backg.mat');
tmp2 = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/dataset_error_with_backg.mat');
wni = 960; wnf = 1300;
idx1 = find(findClosest(wni, tmp1.wavenumber));
idx2 = find(findClosest(wnf, tmp1.wavenumber));
wn = tmp1.wavenumber(idx1:idx2);
Mdatos = tmp1.X(2:end, idx1:idx2)';
idx3 = find(findClosest(wni, tmp2.wavenumber));
idx4 = find(findClosest(wnf, tmp2.wavenumber));
tmp3 = tmp2.X(:, idx3:idx4);
tmp3 = msbackadj(wn, tmp3', 'SHOWPLOT', 3);
Merror = tmp3 - min(tmp3, [], 1);

% datos del tiempo
tiempo = tmp1.time(2:end);

% x-axis
nu = wn;

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
