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

