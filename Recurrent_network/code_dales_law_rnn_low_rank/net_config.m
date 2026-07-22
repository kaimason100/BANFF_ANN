clearvars
close all
clc
%% Configure the standard recurrent BANFF network used in Figures 1--3.
% The random seed and distributions below match the accepted manuscript:
% N = 6000, eta entries in {-4, 0, 4} with probabilities {0.25, 0.5,
% 0.25}, and Gaussian decoder entries with standard deviation 1/sqrt(N).
rng(1, 'twister')

N = 6000;
k = 3;
nx = 1000;
x = 2*rand(k, nx) - 1;
xv = 2*rand(k, 100) - 1;
c = 2*rand(1, nx) - 1;
cv = 2*rand(1, 100) - 1;
win = 5*(2*rand(N, 1) - 1);

eta = 4*round(2*rand(N, k) - 1);
phi = randn(N, k)/sqrt(N);

scriptDir = fileparts(mfilename('fullpath'));
tab = table2cell(readtable(fullfile(scriptDir, 'dynamics_list.xlsx')));
dynamics = tab(:, 1);
dim = cell2mat(tab(:, 2));
input = cell2mat(tab(:, 3)); %#ok<NASGU>
ic = cell2mat(tab(:, 4:7));
nsys = numel(dynamics);
rescale = cell2mat(tab(:, 8));
clear tab scriptDir

save network_configuration.mat -v7.3

