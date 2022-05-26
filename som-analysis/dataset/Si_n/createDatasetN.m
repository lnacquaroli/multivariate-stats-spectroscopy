clc
clear
close all

% Folder
folder_name = '/home/leniac/Work/2019/modos_sup/dataset/Si_n/';
% Read files' names
f = dir( [folder_name, '*.dat'] );
tmp1 = load([folder_name, f(1).name]);
s.lambda = tmp1(:,1);
data = tmp1(:,2);
for idx = 2:length(f)
    tmp2 = load([folder_name, f(idx).name]);
    data = [data tmp2(:, 2)];
end
s.Xerror = data(:, end-11:end);
s.X = data(:, 1:end-12);
s.time = [0 3 9 21 21 51 111 231 231 351 351 471 471 1011 1971 3411 4851 6291 7731];

save('datasetSi_n.mat', '-struct', 's')