function [Results] = target_transformation_pca(spectra_target, initial_conc, X, S, k)
% target_transformation_pca.m
% To be used with factor_analysis_pca.m to determine the real component spectra
%
% Input:
%    spectra_target = Spectra for which you have the targets, Vector (e.g., [1; 10])
%    initial_conc = Initial concentration of the selected spectra, Vector (e.g., [0; 1])
%    X = Matrix with data
%    S = Structure as results from factor_analysis_pca.m
%    k = number of factors
% 
%
% Output:
%
% Results (structure), with fields:
%   R = Row matrix
%   D = Diagonal matrix
%   T = Transformation matrix
%   RFM = Real factors matrix
%   RCM = Real concentration matrix
%   Xbar = Reproduced real spectra
%   rms = rms error between s.dr and Xbar
%
% author: leniac
%

% Selected spectra
spectra_target = spectra_target(:);
spectra = X(:, spectra_target);

% Estimate final concentration
initial_conc = initial_conc(:);
final_conc = 1 - sum(initial_conc, 2);

% Concentration matrix
conc_matrix = [initial_conc final_conc];
conc_matrix = pinv(conc_matrix);

% Row matrix
Row_matrix = S.u * S.s;

% Diagonal matrix
Diag_matrix = diag(S.ev(1:k).^(-1));

% Transformation matrix
tmp1 = spectra * conc_matrix;
Transformation_matrix = Diag_matrix * Row_matrix' * tmp1;

% Real factor matrix
Real_factors_matrix = Row_matrix * Transformation_matrix;

% Real concentration matrix
Real_conc_matrix = Transformation_matrix^(-1) * S.v';

% Real spectra reproduction
Real_spectra_reprod = Real_factors_matrix * Real_conc_matrix;

% Error in reproduction
rms = norm(S.dr - Real_spectra_reprod) / norm(S.dr);

% Wrap up results
Results = struct('R', Row_matrix,...
                 'D', Diag_matrix,...
                 'T', Transformation_matrix,...
                 'RFM', Real_factors_matrix,...
                 'RCM', Real_conc_matrix,...
                 'Xbar', Real_spectra_reprod,...
                 'rms', rms);

end