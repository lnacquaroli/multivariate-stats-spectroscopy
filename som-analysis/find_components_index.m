function [ind_comp, labels_graphs] = find_components_index(comp_names)

wn_comp = [450 626 668 844 874 912 946 986 1060 1146]';
% ind_comp = [75 120 158 180 194 213 228 318 338 424];

ind_comp = zeros(size(wn_comp, 1), 1);
for i = 1:length(ind_comp)
    tmp1 = find(strcmp(num2str(wn_comp(i)), comp_names));
    k = 1;
    if isempty(tmp1)
        while isempty(tmp1)
            tmp1 = find(strcmp(num2str(wn_comp(i)+k), comp_names));
            k = k + 1;
        end
    end
    ind_comp(i) = tmp1;
end

labels_graphs = comp_names(ind_comp);

end


