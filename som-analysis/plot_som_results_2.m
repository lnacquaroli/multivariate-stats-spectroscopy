function [hits] = plot_som_results(mode, Xmaps, Xsom, map2label, ind_comp, cmap)
%plot_som_results Summary of this function goes here

% Add labels to selected (best) map
Xmaps(map2label) = som_label(Xmaps(map2label), 'clear', 'all');
Xmaps(map2label) = som_autolabel(Xmaps(map2label), Xsom, 'all');
sMap = Xmaps(map2label);

switch mode
    case 'umat'
        % U-matrix
        figure()
        som_show(sMap, 'umat', 'all', 'footnote', '', 'colormap', cmap);
        xlabel(['U-matrix, #' num2str(map2label)])
        hits = [];
    case 'bmu'
        % BMU
        plot_bmu(sMap, map2label)
        hits = [];
    case 'hits'
        % Hits
        c = colors_palettes();
        c = c.pure_argyleelegance;
        figure()
        hits = som_hits(sMap, Xsom, 'crisp');
        som_show(sMap, 'empty', 'Hits', 'footnote', '');
        som_show_add('hit', hits, 'Marker', 'lattice', 'MarkerColor', c(2,:), 'Text', 'on',...
             'TextColor', [0 0 0], 'TextSize', 11);
        xlabel(['Hits, #' num2str(map2label)])
    case 'comp_plane'
        % Component planes
        figure()
        som_show(sMap, 'umat', 'all', 'comp', ind_comp, 'norm', 'd', 'footnote', '', 'colormap', cmap);
        hits = [];
    case 'comp_plane_traj'
        % Component planes and trajectory and bmus
        c = my_palette();
        bmus = som_bmus(sMap, Xsom);
        figure()
        h = som_show(sMap, 'umat', 'all', 'comp', ind_comp, 'norm', 'd', 'footnote', '', 'colormap', cmap);
        for i=1:length(ind_comp)
            som_show_add('traj', bmus, 'Markercolor', c(2,:), 'MarkerSize', 5,...
                             'TrajColor', c(2,:), 'TrajWidth', 1.0, 'WidthFactor', 'hit', ...
                             'SizeFactor', 'hit', 'Subplot', i+1);
            som_show_add('label', sMap, 'Textsize', 8, 'TextColor', c(1,:), 'Subplot', i+1);
        end
        % xlabel(['Map, #' num2str(map2label)])
        h.label(1).String = ['U-matrix, Map #' num2str(map2label)];
        hits = [];
end
end

