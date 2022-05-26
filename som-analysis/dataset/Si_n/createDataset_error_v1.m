% Create dataset for the Si-n
clc
clear
close all

% Folder
folder_name = '/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/';
% Read files' names
f = dir( [folder_name, '*.dat'] );
tmp1 = load([folder_name, f(1).name]);
data = tmp1;
for idx = 2:length(f)
    tmp2 = load([folder_name, f(idx).name]);
    data = [data tmp2(:, 2)];
end
data(:, 2:21) = [];
data = sortrows(data, 1);

figure()
plot(data(:,1), data(:,2)), hold on
plot(data(:,1), data(:,end))

% Remove negative values
% tmp0(:,2:end) = tmp0(:,2:end) + abs(min(tmp0(:,2:end),[],2));

% Background correction
% tmp3 = msbackadj(data(:,1), data(:, 2:end), 'SHOWPLOT', 3);
% tmp3 = msbackadj(tmp0(:,1), tmp3, 'SHOWPLOT', 3);

s.wavenumber = data(:,1);
s.X = data(:, 2:end)';
s.X(s.X<0) = 0;
% s.time = [0 3 9 21 51 111 231 471 1011 1971 3411 4851 6291];
s.data = [s.wavenumber, s.X'];
% s.labels = cell(split(num2str(s.time)));
s.comp_names = cell(split(num2str(s.wavenumber')));

figure()
plot(s.wavenumber, s.X(1,:)), hold on
plot(s.wavenumber, s.X(end,:))

save('dataset_error_with_backg.mat', '-struct', 's')
% M = s.X;
% save dataset_error_with_backg.dat M -ascii
