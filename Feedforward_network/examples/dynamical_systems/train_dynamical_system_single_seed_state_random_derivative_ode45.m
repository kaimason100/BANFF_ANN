%% Train one seeded derivative-field dynamical-system network locally.
% This single-network wrapper calls the same shared implementation used by
% the seeded batch trainer. Use it for targeted reruns with a chosen task, seed,
% learning rate, epoch count, and output network set.

clear; close all; clc;

repoRoot = locateRepoRoot();
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'examples')));

set(0, 'DefaultFigureVisible', 'on');

TASK_NAME = "Lorenz";
SEED = 0;
TOTAL_EPOCHS = 100000;
LEARNING_RATE = 0.01;
OUTPUT_NETWORK_SET = sprintf('seeded_state_random_derivative_ode45_lr_%s', learningRateLabel(LEARNING_RATE));

taskRows = containers.Map( ...
    {'Lorenz', 'MO0', 'MO5', 'MO7', 'MO13', 'Rikitake', 'SprottB', 'SprottC', 'SprottS'}, ...
    {5, 34, 39, 41, 47, 32, 8, 9, 25});

taskKey = char(TASK_NAME);
assert(isKey(taskRows, taskKey), 'Unknown dynamical-system task: %s.', taskKey);

trainingOptions.TotalEpochs = TOTAL_EPOCHS;
trainingOptions.ValidationEveryEpochs = 200;
trainingOptions.ValidationSamples = 100;
trainingOptions.ValidationRandomSeed = 9000;
trainingOptions.ScalingSampleTime = 1000;
trainingOptions.RandomSamplesPerEpoch = 1000;
trainingOptions.SampleRandomSeed = 12345;
trainingOptions.StateRandomRange = [-1 1];
trainingOptions.SaveProgressFigures = true;
trainingOptions.ProgressEveryEpochs = 100;
trainingOptions.LearnRate = LEARNING_RATE;
trainingOptions.AdamGradientDecayFactor = 0.99;
trainingOptions.AdamSquaredGradientDecayFactor = 0.999;
trainingOptions.AdamEpsilon = 5e-6;
trainingOptions.OutputNetworkSet = OUTPUT_NETWORK_SET;

fprintf('[local derivative-field single] task=%s | seed=%d | epochs=%d | learning rate=%.6g | output set=%s\n', ...
    taskKey, double(SEED), double(TOTAL_EPOCHS), double(LEARNING_RATE), OUTPUT_NETWORK_SET);
run_dynamical_system_state_random_derivative_ode45(taskRows(taskKey), double(SEED), trainingOptions);

function label = learningRateLabel(value)
label = strrep(sprintf('%.15g', double(value)), '.', 'p');
label = strrep(label, '-', 'm');
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
