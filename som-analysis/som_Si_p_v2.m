% Plot data from som_Si_p_v0.m

clc
clear
% close all

% Choose colormap
% cmap = colormap(brighten(cmocean('amp'), 0.15)); %, 'pivot', 1.43;
% close()
% [brmap, ~, ~] = brewermap(256, 'YlOrRd');
[brmap, ~, ~] = brewermap(256, 'Greys');
cmap = brighten(brmap, 0.15);


% Load data
% load('workspace_som_Si_p_norm_var_v0.mat', 'Xmaps', 'fit_parameters', 'Xsom', 'X', 'nonempty_col')
load('workspace_som_Si_p_norm_lin_nomsize_v0.mat', 'Xmaps', 'fit_parameters', 'Xsom', 'X', 'nonempty_col')

% Identify peaks
figure(), hold on
plot(X.wavenumber, X.X(1,:))
plot(X.wavenumber, X.X(end,:))

% Get indices and labels of components of interest
[ind_comp, labels_graphs] = find_components_index(Xsom.comp_names);

% Best map
% [~, map2label] = min(fit_parameters.MQE, [], 1); % carefull if there are more than one
tmp1 = sortrows(fit_parameters,2);
map2label = table2array(tmp1(2, 1));

%% Make plots
% U-matrix
[~] = plot_som_results('umat', Xmaps, Xsom, map2label, [], cmap);
% BMU
[~] = plot_som_results('bmu', Xmaps, Xsom, map2label, [], []);
% Hits
hits = plot_som_results('hits', Xmaps, Xsom, map2label, [], []);
% CP
[~] = plot_som_results('comp_plane', Xmaps, Xsom, map2label, ind_comp, cmap);
% CP and trajectories
[~] = plot_som_results('comp_plane_traj', Xmaps, Xsom, map2label, ind_comp, cmap);
