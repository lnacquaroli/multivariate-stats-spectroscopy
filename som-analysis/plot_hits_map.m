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

