function [Xmap, fit_parameters] = run_som_maps(Xsom, num_lininit_maps, num_rand_maps, msize)


if isempty(msize) % Default choice
    for i = 1:num_lininit_maps
        Xmap(i) = som_make(Xsom);%, 'msize', [8 5]);
    end
else
    for i = 1:num_lininit_maps
        Xmap(i) = som_make(Xsom, 'msize', msize);
    end
end

mapsize = 1; % lininit map
for i = num_lininit_maps+1:num_lininit_maps+num_rand_maps
    Xmap(i) = som_make(Xsom, 'init', 'randinit', 'seq', 'msize', Xmap(mapsize).topol.msize,...
                       'training', [100 5000]);
end

n = length(Xmap);
fit_prms = zeros(n, 3);
for i = 1:n;
    fit_prms(i,1) = i;
    [fit_prms(i,2), fit_prms(i,3)] = som_quality(Xmap(i), Xsom);
end
fit_parameters = array2table(fit_prms);
fit_parameters.Properties.VariableNames = {'Run', 'MQE', 'TGE'}

end

