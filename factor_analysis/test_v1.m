%% author: leniac
close all
clear
clc

%% Load data

tmp1 = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/datasetSi_n_with_backg.mat');
tmp2 = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/dataset_error_with_backg.mat');

% Select wavenumber range 960-1300 cm^{-1}
wni = 960; wnf = 1300;

% Data array
idx1 = find(findClosest(wni, tmp1.wavenumber));
idx2 = find(findClosest(wnf, tmp1.wavenumber));
wn = tmp1.wavenumber(idx1:idx2);
X = tmp1.X(2:end, idx1:idx2)';
figure, plot(wn,X(:,[1 end]))

% Error array
idx3 = find(findClosest(wni, tmp2.wavenumber));
idx4 = find(findClosest(wnf, tmp2.wavenumber));
tmp3 = tmp2.X(:, idx3:idx4);
% Background correction
tmp3 = msbackadj(wn, tmp3', 'SHOWPLOT', 3);
% Non-negative values
Xerror = tmp3 - min(tmp3, [], 1);
figure, plot(wn,Xerror)

% Samples variable (time, etc)
x1 = tmp1.time(2:end);
%x1 = x1 - x1(1);

% Variation of a spectrum (wavelength, wavenumber, etc)
x2 = wn;

%% Principal factor analysis preliminary: determine the number of components
% If you know them beforehand ignore this section and select k

% PFA to the error matrix
[~] = factor_analysis_pca(Xerror);

% PFA to the data matrix
[~] = factor_analysis_pca(X);

% Number of fators chosen
k = 2;

%% PFA to the data matrix with k factors: abstract 
results1 = factor_analysis_pca(X, k);

%% Target transformation: real 

% Spectra for which you have the targets
spectra_target = [1; size(X, 2)];

% Initial concentration of the selected spectra
ci = [0.35; 0];

% Perform transformation
results2 = target_transformation_pca(spectra_target, ci, X, results1, k);

%% Plots

wnstr = 'Wavenumber [cm^{-1}]';
abstrs = 'Absorbance [a. u.]';
timestr = 'Time [min]';

figure()
plot(x2, results1.dr)
title('Reproduced abstract data'), ylabel(abstrs), xlabel(wnstr)

figure()
plot(x2, results2.RFM, '-')
title('Real principal components'), ylabel(abstrs), xlabel(wnstr)

figure()
plot(x1, results2.RCM', 'o')
title('Real principal concentrations'), ylabel(abstrs), xlabel(timestr)

figure()
plot(x2, results2.Xbar)
title(sprintf('Reproduced real data, RMS = %.1f%s', results2.rms*100, '%'))
ylabel(abstrs), xlabel(wnstr)


