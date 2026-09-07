%% Add the feedforward package to the MATLAB path
% Run once after opening a fresh MATLAB session. Individual Live Scripts and
% the Slurm launcher also run this setup. No data are loaded and no RNG changes.
feedforwardRoot = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(feedforwardRoot,'src')));
addpath(genpath(fullfile(feedforwardRoot,'examples')));
