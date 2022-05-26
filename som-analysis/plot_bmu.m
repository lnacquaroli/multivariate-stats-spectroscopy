function plot_bmu(sMap, map2label)
% BMU
c = [0 118 255] / 255;
figure()
som_show(sMap, 'empty', '', 'footnote', '');
som_show_add('label', sMap, 'Textsize', 8, 'TextColor', c);
xlabel(['BMU, #' num2str(map2label)])
end

