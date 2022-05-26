%% author: leniac

close all
clear
clc

%% Load data

% Input error matrix: usually the same sample measured at least 10 times consecutively
Merror = load('MatrizError.dat');

% Input data matrix
Mdatos = load('MatrizDatos.dat');

% Samples variable (time, etc)
x1 = load('ti2.dat');
x1 = x1 - x1(1);

% Variation of a spectrum (wavelength, wavenumber, etc)
x2 = 1:size(Merror,1);

%% Principal factor analysis preliminary: determine the number of components
% If you know them beforehand ignore this section and select k

% PFA to the error matrix
[~] = factor_analysis_pca(Merror);

% PFA to the data matrix
[~] = factor_analysis_pca(Mdatos);

% Number of fators chosen
k = 2;

%% PFA to the data matrix with k factors: abstract 
results1 = factor_analysis_pca(Mdatos, k);

%% Target transformation: real 

% Spectra for which you have the targets
spectra_target = [1; size(Mdatos, 2)];

% Initial concentration of the selected spectra
ci = [1; 0];

% Perform transformation
results2 = target_transformation_pca(spectra_target, ci, Mdatos, results1, k);

%% Plots

figure()
plot(x2, results1.dr)
title('Reproduced abstract data'), ylabel('Data [units]'), xlabel('x2 [units]')

figure()
plot(x2, results2.RFM, '-')
title('Real principal components'), ylabel('[units]'), xlabel('x2 [units]')

figure()
plot(x1, results2.RCM', 'o')
title('Real principal concentrations'), ylabel('[units]'), xlabel('x1 [units]')

figure()
plot(x2, results2.Xbar)
title(sprintf('Reproduced real data, RMS = %.1f%s', results2.rms*100, '%')), ylabel('[units]'), xlabel('x2 [units]')


