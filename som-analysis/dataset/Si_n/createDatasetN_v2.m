% Create dataset for the Si-n
clc
clear
close all

% Folderdfd
tmp0 = load('/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/all_spectra.txt');
% folder_name = '/home/leniac/Work/2019/modos_sup/som_analysis/dataset/Si_n/';
% Read files' names
% f = dir( [folder_name, '*.dat'] );
% tmp1 = load([folder_name, f(1).name]);
% data = tmp1;
% for idx = 2:length(f)
%     tmp2 = load([folder_name, f(idx).name]);
%     data = [data tmp2(:, 2)];
% end

tmp0(:, 9) = [];
figure()
plot(tmp0(:,1), tmp0(:,2)), hold on
plot(tmp0(:,1), tmp0(:,end))

% Background correction
tmp3 = msbackadj(tmp0(:,1), tmp0(:, 2:end), 'SHOWPLOT', 3);
% tmp3 = msbackadj(tmp0(:,1), tmp3, 'SHOWPLOT', 3);

s.wavenumber = tmp0(:,1);
s.X = tmp3';
s.time = [0 3 9 21 51 111 231 471 1011 1971 3411 4851 6291];
s.data = [s.wavenumber, s.X'];
s.labels = cell(split(num2str(s.time)));
s.comp_names = cell(split(num2str(s.wavenumber')));
% Remove negative values
s.X = s.X + abs(min(s.X,[],2));

figure()
plot(s.wavenumber, s.X(1,:)), hold on
plot(s.wavenumber, s.X(end,:))

save('datasetSi_n.mat', '-struct', 's')
