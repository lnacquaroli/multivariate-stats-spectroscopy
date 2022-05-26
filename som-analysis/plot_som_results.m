function [hits] = plot_som_results(mode, Xmaps, Xsom, map2label, ind_comp, cmap)
%plot_som_results Summary of this function goes here

% Add labels to selected (best) map
sMap = map_autolabel(Xmaps, Xsom, map2label);

switch mode
    case 'umat'
        % U-matrix
        plot_umatrix(sMap, map2label, cmap)
        hits = [];
    case 'bmu'
        % BMU
        plot_bmu(sMap, map2label)
        hits = [];
    case 'hits'
        % Hits
        hits = plot_hits_map(sMap, Xsom, map2label);
    case 'comp_plane'
        figure()
        % Component planes
        plot_cplanes(sMap, ind_comp, cmap)
        hits = [];
    case 'comp_plane_traj'
        % Component planes and trajectory and bmus
        plot_cp_trajectory_bmu(sMap, Xsom, ind_comp, cmap, map2label)
        hits = [];
end
end

function [sMap] = map_autolabel(Xmaps, Xsom, map2label)
% Add labels to map
Xmaps(map2label) = som_label(Xmaps(map2label), 'clear', 'all');
Xmaps(map2label) = som_autolabel(Xmaps(map2label), Xsom, 'all');
sMap = Xmaps(map2label);
end

function plot_umatrix(sMap, map2label, cmap)
% U-matrix
figure()
som_show(sMap, 'umat', 'all', 'footnote', '', 'colormap', cmap);
xlabel(['U-matrix, #' num2str(map2label)])
end

function plot_bmu(sMap, map2label)
% BMU
c = [0 118 255] / 255;
figure()
som_show(sMap, 'empty', '', 'footnote', '');
som_show_add('label', sMap, 'Textsize', 8, 'TextColor', c);
xlabel(['BMU, #' num2str(map2label)])
end

function [hits] = plot_hits_map(sMap, Xsom, map2label)
% BMU
c = colors_palettes();
c = c.pure_argyleelegance;
figure()
hits = som_hits(sMap, Xsom, 'crisp');
som_show(sMap, 'empty', 'Hits', 'footnote', '');
som_show_add('hit', hits, 'Marker', 'lattice', 'MarkerColor', c(2,:), 'Text', 'on',...
             'TextColor', [0 0 0], 'TextSize', 11);
xlabel(['Hits, #' num2str(map2label)])
end

function plot_cp_trajectory_bmu(sMap, Xsom, ind_comp, cmap, map2label)
c = my_palette();
bmus = som_bmus(sMap, Xsom);
figure();
h = som_show(sMap, 'umat', 'all', 'comp', ind_comp, 'norm', 'd', 'footnote', '', 'colormap', cmap);
for i=1:length(ind_comp)
    som_show_add('traj', bmus, 'Markercolor', c(2,:), 'MarkerSize', 5,...
                 'TrajColor', c(2,:), 'TrajWidth', 1.0, 'WidthFactor', 'hit', ...
                 'SizeFactor', 'hit', 'Subplot', i+1);
    som_show_add('label', sMap, 'Textsize', 8, 'TextColor', c(1,:), 'Subplot', i+1);
end
% xlabel(['Map, #' num2str(map2label)])
h.label(1).String = ['U-matrix, Map #' num2str(map2label)];
end

function plot_cplanes(sMap, ind_comp, cmap)
figure()
som_show(sMap, 'umat', 'all', 'comp', ind_comp, 'norm', 'd', 'footnote', '', 'colormap', cmap);
end
