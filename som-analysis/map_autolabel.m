function [sMap] = map_autolabel(Xmaps, Xsom, map2label)
%map_autolabel
Xmaps(map2label) = som_label(Xmaps(map2label), 'clear', 'all');
Xmaps(map2label) = som_autolabel(Xmaps(map2label), Xsom, 'all');
sMap = Xmaps(map2label);
end

