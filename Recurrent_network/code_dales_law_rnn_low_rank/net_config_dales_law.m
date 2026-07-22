clearvars
close all
clc
%% Configure the Dale-constrained recurrent BANFF network used in Figure S4.
% The network has 3000 excitatory and 3000 inhibitory neurons. Encoder
% entries are 0 or 4 with equal probability. Decoder rows carry the sign of
% the presynaptic neuron and use rounded Gaussian components scaled by
% 1/sqrt(N), matching the accepted manuscript.
rng(1, 'twister')

N = 6000;
k = 3;
nx = 1000;
x = 2*rand(k, nx) - 1;
xv = 2*rand(k, 100) - 1;
c = 2*rand(1, nx) - 1;
cv = 2*rand(1, 100) - 1;
win = 5*(2*rand(N, 1) - 1);

NE = N/2;
eta = 4*round(rand(N, k));
phiMagnitude = abs(round(randn(N, k)))/sqrt(N);
phi = phiMagnitude;
phi((NE + 1):end, :) = -phiMagnitude((NE + 1):end, :);

scriptDir = fileparts(mfilename('fullpath'));
tab = table2cell(readtable(fullfile(scriptDir, 'dynamics_list.xlsx')));
dynamics = tab(:, 1);
dim = cell2mat(tab(:, 2));
input = cell2mat(tab(:, 3)); %#ok<NASGU>
ic = cell2mat(tab(:, 4:7));
nsys = numel(dynamics);
rescale = cell2mat(tab(:, 8));
clear tab scriptDir phiMagnitude

save network_configuration_dales_law.mat -v7.3
