%% Continue selected seeded derivative-field dynamical-system networks locally.
% This wrapper uses the same shared trainer as the seeded batch implementation.
% It loads selected networks from the source network set, continues bias-only
% derivative-field training, and saves into a separate continuation network set.

clear; close all; clc;

repoRoot = locateRepoRoot();
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'examples')));

set(0, 'DefaultFigureVisible', 'on');

% Edit this list to choose the exact networks to continue locally.
% Rows are {taskName, seed}.
continuationPairs = {
    'Lorenz', 5
    };

% Extra epochs for the continuation stage.
CONTINUATION_EPOCHS = 50000;

taskRowsByName = containers.Map( ...
    {'Lorenz', 'MO0', 'MO5', 'MO7', 'MO13', 'Rikitake', 'SprottB', 'SprottC', 'SprottS'}, ...
    {5, 34, 39, 41, 47, 32, 8, 9, 25});

trainingOptions.TotalEpochs = CONTINUATION_EPOCHS;
trainingOptions.ValidationEveryEpochs = 200;
trainingOptions.ValidationSamples = 100;
trainingOptions.ValidationRandomSeed = 9000;
trainingOptions.ScalingSampleTime = 1000;
trainingOptions.RandomSamplesPerEpoch = 1000;
trainingOptions.SampleRandomSeed = 12345;
trainingOptions.StateRandomRange = [-1 1];
trainingOptions.SaveProgressFigures = true;
trainingOptions.ProgressEveryEpochs = 100;
trainingOptions.LearnRate = 0.001;
trainingOptions.AdamGradientDecayFactor = 0.99;
trainingOptions.AdamSquaredGradientDecayFactor = 0.999;
trainingOptions.AdamEpsilon = 5e-6;
trainingOptions.ContinueFromSource = true;
trainingOptions.SourceNetworkSet = 'seeded_state_random_derivative_ode45_lr_0p01';
trainingOptions.OutputNetworkSet = sprintf('seeded_state_random_derivative_ode45_continued_lr_0p001_extra_%d', double(trainingOptions.TotalEpochs));
trainingOptions.ContinuationLabel = sprintf('continued_lr_0p001_extra_%d', double(trainingOptions.TotalEpochs));

fprintf('[local derivative-field continuation] source set=%s | output set=%s\n', trainingOptions.SourceNetworkSet, trainingOptions.OutputNetworkSet);
fprintf('[local derivative-field continuation] extra epochs=%d | learning rate=%.6g | validation every %d epochs | validation samples=%d\n', ...
    trainingOptions.TotalEpochs, trainingOptions.LearnRate, trainingOptions.ValidationEveryEpochs, trainingOptions.ValidationSamples);
fprintf('[local derivative-field continuation] training samples=%d | sample RNG seed=%d | validation RNG seed=%d\n', ...
    trainingOptions.RandomSamplesPerEpoch, trainingOptions.SampleRandomSeed, trainingOptions.ValidationRandomSeed);

for pairIdx = 1:size(continuationPairs, 1)
    taskName = char(continuationPairs{pairIdx, 1});
    weightSeed = double(continuationPairs{pairIdx, 2});
    assert(isKey(taskRowsByName, taskName), 'Unknown continuation task: %s.', taskName);
    assert(isscalar(weightSeed) && isfinite(weightSeed) && weightSeed == round(weightSeed), ...
        'Continuation seed must be a scalar integer for %s.', taskName);
    fprintf('\n[local derivative-field continuation] pair %d/%d | task=%s | seed=%d\n', ...
        pairIdx, size(continuationPairs, 1), taskName, weightSeed);
    run_dynamical_system_state_random_derivative_ode45(taskRowsByName(taskName), weightSeed, trainingOptions);
end

function repoRoot = locateRepoRoot()
starts = string({pwd, fileparts(mfilename('fullpath'))});
for s = starts
    candidate = char(s);
    while ~isempty(candidate)
        if isfolder(fullfile(candidate, 'trained_networks')) && isfolder(fullfile(candidate, 'examples'))
            repoRoot = candidate;
            return
        end
        parent = fileparts(candidate);
        if strcmp(parent, candidate), break; end
        candidate = parent;
    end
end
error('Could not locate repository root.');
end
