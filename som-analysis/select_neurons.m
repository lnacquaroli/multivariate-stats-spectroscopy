function [REF_eem, neurons_int, neurons_int_non_sort] = select_neurons(som_map, data_raw, hit1, n)
%select_neurons Summary of this function goes here

%% View neurons of interest
% Retrieve the spectral information stored in the neuron of interest
EEM_den = som_denormalize(som_map);
% EEM_den = som_map;
[~, I] = sort(hit1, 'descend'); % with 'var' normalization
select_neurons = 1:n;
neurons_int = sort(I(select_neurons));
% REF_vec = zeros(length(neurons_int), length(nonempty));
REF_vec = zeros(length(neurons_int), length(data_raw.Em)*length(data_raw.Ex));
REF_eem = zeros(length(neurons_int), length(data_raw.Em), length(data_raw.Ex));
figure(), xlabel('Ex. (nm)'), ylabel('Em. (nm)')
for i = 1:length(neurons_int)
    REF_vec(i, :) = EEM_den.codebook(neurons_int(i), :);
    REF_eem(i, :, :) = reshape(REF_vec(i, :), length(data_raw.Em), length(data_raw.Ex));
    subplot(3,4,i), contourf(data_raw.Ex, data_raw.Em, squeeze(REF_eem(i, :, :))), colorbar
    title(['Neuron ' num2str(neurons_int(i))])
end

neurons_int_non_sort = I(select_neurons);

end

