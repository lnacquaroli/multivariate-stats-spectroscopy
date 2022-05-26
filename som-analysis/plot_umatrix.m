function plot_umatrix(sMap, map2label, cmap)
% U-matrix
figure()
som_show(sMap, 'umat', 'all', 'footnote', '', 'colormap', cmap);
xlabel(['U-matrix, #' num2str(map2label)])
end

