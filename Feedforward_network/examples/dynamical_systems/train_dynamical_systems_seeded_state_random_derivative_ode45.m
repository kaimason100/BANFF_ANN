%% Train seeded derivative-field dynamical-system networks locally.
% This wrapper calls the shared derivative-field trainer with the publication
% default options. For a given task/seed, the saved network, metadata, training
% samples, validation samples, and output folder are deterministic up to normal
% MATLAB numerical differences across machines.

clear; close all; clc;

repoRoot = locateRepoRoot();
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'examples')));

set(0, 'DefaultFigureVisible', 'on');

% Full local batch. Edit these to run a subset.
SELECTED_TASKS = ["Lorenz", "MO0", "MO5", "MO7", "MO13", "Rikitake", "SprottB", "SprottC", "SprottS"];
SEEDS = 0:9;

% These defaults intentionally match local derivative-field batch implementation.
trainingOptions.TotalEpochs = 100000;
trainingOptions.ValidationEveryEpochs = 200;
trainingOptions.ValidationSamples = 100;
trainingOptions.ValidationRandomSeed = 9000;
trainingOptions.ScalingSampleTime = 1000;
trainingOptions.RandomSamplesPerEpoch = 1000;
trainingOptions.SampleRandomSeed = 12345;
trainingOptions.StateRandomRange = [-1 1];
trainingOptions.SaveProgressFigures = true;
trainingOptions.ProgressEveryEpochs = 100;

% Adam optimiser parameters.
trainingOptions.LearnRate = 0.01;
trainingOptions.AdamGradientDecayFactor = 0.99;
trainingOptions.AdamSquaredGradientDecayFactor = 0.999;
trainingOptions.AdamEpsilon = 5e-6;

taskRows = containers.Map( ...
    {'Lorenz', 'MO0', 'MO5', 'MO7', 'MO13', 'Rikitake', 'SprottB', 'SprottC', 'SprottS'}, ...
    {5, 34, 39, 41, 47, 32, 8, 9, 25});

for taskName = SELECTED_TASKS
    taskKey = char(taskName);
    assert(isKey(taskRows, taskKey), 'Unknown dynamical-system task: %s.', taskKey);
    for seedValue = SEEDS
        fprintf('\n[local derivative-field] Running %s seed %d\n', taskKey, double(seedValue));
        run_dynamical_system_state_random_derivative_ode45(taskRows(taskKey), double(seedValue), trainingOptions);
    end
end

function repoRoot = locateRepoRoot()
starts = string({pwd, fileparts(mfilename('fullpath'))});
for s = starts
    candidate = char(s);
    while ~isempty(candidate)
        if isfolder(fullfile(candidate, 'examples')) && isfolder(fullfile(candidate, 'src'))
            repoRoot = candidate;
            return
        end
        nestedCandidate = fullfile(candidate, 'Feedforward_network');
        if isfolder(fullfile(nestedCandidate, 'examples')) && isfolder(fullfile(nestedCandidate, 'src'))
            repoRoot = nestedCandidate;
            return
        end
        parent = fileparts(candidate);
        if strcmp(parent, candidate), break; end
        candidate = parent;
    end
end
error('Could not locate repository root.');
end
