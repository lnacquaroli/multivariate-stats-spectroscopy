% Run SOM

clc
clear
close all

% Load data
X = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_p/datasetSi_p.mat');

% Create SOM data struct
Xsom = som_data_struct(X.X, 'labels', X.labels, 'comp_names', X.comp_names);

% Scale data to unit variance
Xsom = som_normalize(Xsom, 'range');
% Xsom = som_normalize(Xsom, 'var');
% plot(X.wavenumber, X.X(1,:)), hold on
% plot(X.wavenumber, X.X(6,:))
% plot(X.wavenumber, X.X(end,:))

% Remove zero-columns
empty_col = find(sum(Xsom.data) == 0);
nonempty_col = find(sum(Xsom.data)~= 0); 
Xsom.data(:, empty_col) = [];
Xsom.comp_names(empty_col) = [];
Xsom.comp_norm(empty_col) = [];

% Run SOM for several times
[Xmaps, fit_parameters] = run_som_maps(Xsom, 1, 20, [7 3]);

save workspace_som_Si_p_norm_lin_nomsize_v0.mat % same as Si-n in fact
% save workspace_som_Si_p_norm_var_v0.mat
