% Create dataset for the Si-p 

clc
clear
close all

% Folder
folder_name = '/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_p/';
% Read files' names
f = dir( [folder_name, '*.txt'] );
tmp1 = load([folder_name, f(1).name]);
s.wavenumber = tmp1(:,1);
data = tmp1(:,2);
for idx = 2:length(f)
    tmp2 = load([folder_name, f(idx).name]);
    data = [data tmp2(:, 2)];
end
% s.Xerror = data(:, end-11:end);
s.X = data(:,1:end-1)';
s.time = [0 3 8 15 30 60 120 235 465 970 1498]; % min

% Remove negative values
s.X = s.X + abs(min(s.X,[],2));

% Sort data for background correction
s.data = [s.wavenumber, s.X'];
s.data = sortrows(s.data, 1);

% Background correction
% tmp3 = msbackadj(s.data(:,1), s.data(:, 2:end), 'SHOWPLOT', 3);
% tmp3 = msbackadj(s.data(:,1), tmp3, 'SHOWPLOT', 3);

% Save data and zero non-negatives
s.X = s.data(:, 2:end)';
% s.X(s.X<0) = 0;
s.data = [s.data(:,1), s.X'];
s.wavenumber = s.data(:,1);
s.labels = cell(split(num2str(s.time)));
s.comp_names = cell(split(num2str(s.wavenumber')));

plot(s.wavenumber, s.X(1,:)), hold on
plot(s.wavenumber, s.X(end,:))

save('datasetSi_p_with_backg.mat', '-struct', 's')

