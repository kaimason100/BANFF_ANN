%% Test seeded state-random derivative-field dynamical-system networks.
% This script loads networks trained by run_dynamical_system_state_random_derivative_ode45
% and tests them by integrating the learned continuous-time vector field:
%   dx_norm/dt = network(x_norm) - x_norm
% with ode45, matching the reference derivative-field code's testing route.
% Training selects parameters by fixed random derivative-field validation
% loss; this script performs the final closed-loop phase-portrait sliced
% Wasserstein test.

clear; close all; clc;

NETWORK_SET = "seeded_state_random_derivative_ode45_lr_0p01";

% Set this to true to automatically test a continuation-trained network when
% one exists for the current task/seed. Set it to false to force testing of the
% base NETWORK_SET only, even if continuation networks are present.
USE_CONTINUATION_IF_AVAILABLE = true;

% If continuation use is enabled, leave this as "auto" to scan all continuation
% folders and prefer the largest extra-epoch continuation. Or set a string array
% of network-set names to control the search order explicitly.
CONTINUATION_NETWORK_SETS = "auto";
TASKS = ["Lorenz", "MO0", "MO5", "MO7", "MO13", "Rikitake", "SprottB", "SprottC", "SprottS"];
SEEDS = 0:9;

CLOSED_LOOP_ROLLOUT_LENGTH = 1000;
ODE_OUTPUT_DT = 0.01; % ode45 output grid only, matching the reference code; not a fixed-step integrator.
TEST_IC_RANDOM_SEED = 9100;
TEST_IC_PERTURBATION_SCALE = 0; % Standardized-state perturbation around the spreadsheet IC.
PLOT_TEST_FIGURES = true;
SAVE_RESULTS = true;

% The reference code calls ode45 without custom tolerances. Leave these empty
% to use MATLAB's ode45 defaults; tight tolerances can make the large dense
% network RHS prohibitively slow on CPU.
ODE_REL_TOL = [];
ODE_ABS_TOL = [];
ODE_MAX_STEP = [];
PHASE_WASSERSTEIN_OPTIONS = struct( ...
    'NumProjections', 128, ...
    'TrimFraction', 0.10, ...
    'Subsample', 5, ...
    'TransientFraction', 0.10, ...
    'MaxPoints', 1250);

repoRoot = locateRepoRoot();
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(genpath(fullfile(repoRoot, 'examples')));

figureDir = fullfile(repoRoot, 'outputs', 'test_figures', 'dynamical_systems_seeded_state_random_derivative_ode45');
if exist(figureDir, 'dir') ~= 7, mkdir(figureDir); end

results = table('Size', [0, 10], ...
    'VariableTypes', {'string', 'double', 'string', 'logical', 'double', 'double', 'double', 'double', 'double', 'string'}, ...
    'VariableNames', {'Task', 'Seed', 'NetworkSet', 'UsedContinuation', 'PhaseWasserstein', 'BestValidationEpoch', 'BestValidationLoss', 'NumSamples', 'MeanDt', 'NetworkPath'});

for task = TASKS
    seedWD = nan(numel(SEEDS), 1);
    lossHistories = cell(numel(SEEDS), 1);
    for s = 1:numel(SEEDS)
        seedValue = SEEDS(s);
        [modelPath, selectedNetworkSet, usedContinuation] = selectDerivativeModelPath(repoRoot, NETWORK_SET, CONTINUATION_NETWORK_SETS, USE_CONTINUATION_IF_AVAILABLE, task, seedValue);
        assert(isfile(modelPath), 'Missing derivative-field network: %s', modelPath);
        S = load(modelPath);
        assert(isfield(S, 'net') && isa(S.net, 'dlnetwork'), '%s seed %d does not contain a dlnetwork named net.', task, seedValue);
        assert(isfield(S, 'stateStats'), '%s seed %d lacks stateStats.', task, seedValue);
        assert(isfield(S, 'metadata'), '%s seed %d lacks metadata.', task, seedValue);
        assert(isfield(S.metadata, 'trainingSamplingMode') && strcmp(char(S.metadata.trainingSamplingMode), 'state-random-derivative-field'), ...
            '%s seed %d was not trained with the derivative-field implementation.', task, seedValue);
        seedFigureDir = fullfile(figureDir, sprintf('%s_seed_%03d', matlab.lang.makeValidName(char(task)), double(seedValue)));
        if exist(seedFigureDir, 'dir') ~= 7, mkdir(seedFigureDir); end

        relTol = ODE_REL_TOL;
        absTol = ODE_ABS_TOL;
        maxStep = ODE_MAX_STEP;
        tEval = odeOutputTimes(CLOSED_LOOP_ROLLOUT_LENGTH, ODE_OUTPUT_DT);
        x0Norm = chooseStandardizedTestInitialCondition(repoRoot, task, S, TEST_IC_RANDOM_SEED, TEST_IC_PERTURBATION_SCALE);

        odeOptions = makeOdeOptions(relTol, absTol, maxStep);
        fprintf('%s seed %d derivative-field test starting | network set %s | continuation=%d | rollout %.6g | output dt %.6g | evaluator predict\n', ...
            task, seedValue, selectedNetworkSet, usedContinuation, CLOSED_LOOP_ROLLOUT_LENGTH, ODE_OUTPUT_DT);
        testTimer = tic;
        [tTrue, xTrue] = ode45(@(t, x) trueStandardizedDerivative(task, x, S.stateStats), tEval, x0Norm(:), odeOptions);
        [tPred, xPred] = ode45(@(t, x) learnedDerivativeFieldPredict(S.net, x), tEval, x0Norm(:), odeOptions); %#ok<ASGLU>
        elapsedSeconds = toc(testTimer);
        assert(isequal(size(xPred), size(xTrue)), '%s seed %d learned and true trajectories differ in size.', task, seedValue);

        wd = phasePortraitWassersteinDistance(xPred, xTrue, PHASE_WASSERSTEIN_OPTIONS);
        seedWD(s) = wd;
        bestEpoch = metadataNestedField(S, 'validationInfo', 'bestEpoch', NaN);
        bestValLoss = metadataNestedField(S, 'validationInfo', 'bestValidationLoss', NaN);
        meanDt = mean(diff(tTrue), 'omitnan');
        lossHistories{s} = savedLossHistory(S, seedValue);

        fprintf('%s seed %d derivative-field test phase WD %.6g | best validation epoch %.0f | validation loss %.6g | samples %d | mean dt %.6g | elapsed %.2f s\n', ...
            task, seedValue, wd, bestEpoch, bestValLoss, numel(tTrue), meanDt, elapsedSeconds);

        results = [results; {task, seedValue, selectedNetworkSet, usedContinuation, wd, bestEpoch, bestValLoss, numel(tTrue), meanDt, string(modelPath)}]; %#ok<AGROW>

        if PLOT_TEST_FIGURES
            safeName = sprintf('%s_seed_%03d', matlab.lang.makeValidName(char(task)), double(seedValue));
            plotTrajectoryFigure(tTrue, xTrue, xPred, task, seedFigureDir, safeName, wd);
            plotPhaseFigure(xTrue, xPred, task, seedFigureDir, safeName, wd);
            plotReturnMapFigure(xTrue, xPred, task, seedFigureDir, safeName, wd);
        end
        clear S
    end
    fprintf('%s derivative-field phase WD mean %.6g | SD %.6g across %d seeds\n', ...
        task, mean(seedWD, 'omitnan'), std(seedWD, 0, 'omitnan'), sum(isfinite(seedWD)));
    if PLOT_TEST_FIGURES
        plotSystemLossFigure(task, lossHistories, figureDir);
    end
end

disp(results);
perTaskSeedSummary = groupsummary(results, 'Task', {'mean', 'std'}, {'PhaseWasserstein', 'BestValidationLoss'});
disp(perTaskSeedSummary);
distributionFig = figure('Color', 'w');
boxchart(categorical(results.Task), results.PhaseWasserstein);
ylabel('Phase WD'); xlabel('Task'); title('Derivative-field DS phase WD'); grid on;
saveas(distributionFig, fullfile(figureDir, 'dynamical_systems_seeded_state_random_derivative_ode45_wd_distribution.png'));

if SAVE_RESULTS
    outDir = fullfile(repoRoot, 'outputs', 'test_results');
    if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    outFile = fullfile(outDir, 'dynamical_systems_seeded_state_random_derivative_ode45_test_results.mat');
    writetable(results, fullfile(figureDir, 'dynamical_systems_seeded_state_random_derivative_ode45_results.csv'));
    writetable(perTaskSeedSummary, fullfile(figureDir, 'dynamical_systems_seeded_state_random_derivative_ode45_per_task_seed_summary.csv'));
    save(outFile, 'results', 'NETWORK_SET', 'USE_CONTINUATION_IF_AVAILABLE', 'CONTINUATION_NETWORK_SETS', ...
        'TASKS', 'SEEDS', 'CLOSED_LOOP_ROLLOUT_LENGTH', ...
        'ODE_OUTPUT_DT', 'ODE_REL_TOL', 'ODE_ABS_TOL', 'ODE_MAX_STEP', ...
        'TEST_IC_RANDOM_SEED', 'TEST_IC_PERTURBATION_SCALE', 'PHASE_WASSERSTEIN_OPTIONS', 'PLOT_TEST_FIGURES', ...
        'perTaskSeedSummary');
    fprintf('Saved derivative-field test results to %s\n', outFile);
end

function modelPath = seededDerivativeModelPath(repoRoot, networkSet, taskName, seedValue)
folderName = matlab.lang.makeValidName(char(taskName));
modelPath = fullfile(repoRoot, 'trained_networks', char(networkSet), folderName, ...
    sprintf('%s_seed_%03d_network.mat', folderName, double(seedValue)));
end

function [modelPath, selectedNetworkSet, usedContinuation] = selectDerivativeModelPath(repoRoot, baseNetworkSet, continuationNetworkSets, useContinuation, taskName, seedValue)
usedContinuation = false;
selectedNetworkSet = string(baseNetworkSet);
if useContinuation
    candidateSets = resolveContinuationNetworkSets(repoRoot, continuationNetworkSets);
    for k = 1:numel(candidateSets)
        candidatePath = seededDerivativeModelPath(repoRoot, candidateSets(k), taskName, seedValue);
        if isfile(candidatePath)
            modelPath = candidatePath;
            selectedNetworkSet = candidateSets(k);
            usedContinuation = true;
            fprintf('%s seed %d: using continuation network set %s\n', taskName, seedValue, char(selectedNetworkSet));
            return
        end
    end
end
modelPath = seededDerivativeModelPath(repoRoot, baseNetworkSet, taskName, seedValue);
selectedNetworkSet = string(baseNetworkSet);
fprintf('%s seed %d: no continuation network found; using base network set %s\n', taskName, seedValue, char(selectedNetworkSet));
end

function candidateSets = resolveContinuationNetworkSets(repoRoot, continuationNetworkSets)
if isstring(continuationNetworkSets) && isscalar(continuationNetworkSets) && strcmpi(continuationNetworkSets, "auto")
    listing = dir(fullfile(repoRoot, 'trained_networks', 'seeded_state_random_derivative_ode45_continued_lr_0p001_extra_*'));
    listing = listing([listing.isdir]);
    names = string({listing.name});
    if isempty(names)
        candidateSets = strings(0, 1);
        return
    end
    extraEpochs = nan(numel(names), 1);
    for k = 1:numel(names)
        token = regexp(char(names(k)), 'extra_(\d+)$', 'tokens', 'once');
        if ~isempty(token)
            extraEpochs(k) = str2double(token{1});
        end
    end
    folderTimes = [listing.datenum].';
    [~, order] = sortrows([-extraEpochs(:), -folderTimes(:)]);
    candidateSets = names(order);
else
    candidateSets = string(continuationNetworkSets(:));
end
candidateSets = candidateSets(candidateSets ~= "");
end

function x0Norm = chooseStandardizedTestInitialCondition(repoRoot, taskName, S, icSeed, perturbationScale)
dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
row = find(strcmp(string(dyn{:, 1}), string(taskName)), 1, 'first');
assert(~isempty(row), 'Could not find %s in dynamics_list.xlsx.', taskName);
d = dyn{row, 2};
ic = dyn{row, 4:6}.';
baseNorm = applyStateStandardization(ic(1:d).', S.stateStats);
if perturbationScale > 0
    stream = RandStream('mt19937ar', 'Seed', double(icSeed));
    baseNorm = baseNorm + double(perturbationScale) * (2*rand(stream, size(baseNorm)) - 1);
end
x0Norm = baseNorm(:);
end

function value = metadataNestedField(S, structName, fieldName, defaultValue)
value = defaultValue;
if isfield(S, structName) && isfield(S.(structName), fieldName) && ~isempty(S.(structName).(fieldName))
    value = S.(structName).(fieldName);
end
end

function history = savedLossHistory(S, seedValue)
history.seed = double(seedValue);
history.trainingEpoch = [];
history.trainingLoss = [];
history.validationEpoch = [];
history.validationLoss = [];
if isfield(S, 'trainingInfo') && istable(S.trainingInfo)
    T = S.trainingInfo;
    if any(strcmp(T.Properties.VariableNames, 'Epoch')) && any(strcmp(T.Properties.VariableNames, 'TrainingLoss'))
        history.trainingEpoch = double(T.Epoch(:));
        history.trainingLoss = double(T.TrainingLoss(:));
    end
end
if isfield(S, 'validationInfo') && isfield(S.validationInfo, 'history') && istable(S.validationInfo.history)
    V = S.validationInfo.history;
    if any(strcmp(V.Properties.VariableNames, 'Epoch')) && any(strcmp(V.Properties.VariableNames, 'ValidationLoss'))
        history.validationEpoch = double(V.Epoch(:));
        history.validationLoss = double(V.ValidationLoss(:));
    end
end
end

function A = numericOutput(Y)
if isa(Y, 'dlarray')
    A = extractdata(Y);
else
    A = Y;
end
try, A = gather(A); catch, end
A = double(A);
end

function x = applyStateStandardization(x, stats)
x = stats.scale * ((x - stats.mu) ./ stats.sigma);
x(~isfinite(x)) = 0;
end

function xRaw = inverseStateStandardization(xNorm, stats)
xRaw = (double(xNorm) ./ stats.scale) .* stats.sigma + stats.mu;
end

function tEval = odeOutputTimes(rolloutLength, outputDt)
assert(isfinite(rolloutLength) && rolloutLength > 0, 'CLOSED_LOOP_ROLLOUT_LENGTH must be positive.');
assert(isfinite(outputDt) && outputDt > 0, 'ODE_OUTPUT_DT must be positive.');
nWholeSteps = floor(double(rolloutLength) / double(outputDt));
tEval = (0:nWholeSteps).' * double(outputDt);
if tEval(end) < double(rolloutLength)
    tEval(end + 1, 1) = double(rolloutLength);
end
end

function odeOptions = makeOdeOptions(relTol, absTol, maxStep)
odeOptions = odeset();
if ~isempty(relTol)
    odeOptions = odeset(odeOptions, 'RelTol', relTol);
end
if ~isempty(absTol)
    odeOptions = odeset(odeOptions, 'AbsTol', absTol);
end
if ~isempty(maxStep)
    odeOptions = odeset(odeOptions, 'MaxStep', maxStep);
end
end

function dxNormDt = trueStandardizedDerivative(dynamicsName, xNorm, stats)
xNorm = double(xNorm(:));
xRaw = inverseStateStandardization(xNorm.', stats).';
dxRawDt = dynamicsDerivative(dynamicsName, xRaw);
dxNormDt = stats.scale * (dxRawDt(:) ./ stats.sigma(:));
assert(all(isfinite(dxNormDt)), '%s true standardized derivative produced non-finite values.', dynamicsName);
end

function dxdt = learnedDerivativeFieldPredict(net, x)
x = double(x(:));
Y = predict(net, dlarray(single(x), 'CB'));
y = double(numericOutput(Y));
dxdt = y(:) - x;
if any(~isfinite(dxdt)) || any(abs(dxdt) > 1e6)
    dxdt(:) = NaN;
end
end

function dx = dynamicsDerivative(dynamicsName, x)
inputSize = size(x);
if isvector(x)
    x = x(:);
    inputSize = size(x);
end
dx = int_dyn(x, string(dynamicsName), 0, 0, 0, 'Simulate');
assert(isequal(size(dx), inputSize), '%s derivative returned size [%s] for input size [%s].', ...
    dynamicsName, num2str(size(dx)), num2str(inputSize));
end

function distanceValue = phasePortraitWassersteinDistance(xPred, xTrue, options)
if any(~isfinite(xPred(:))) || any(~isfinite(xTrue(:)))
    distanceValue = Inf;
    return
end
nStates = min(size(xPred, 2), size(xTrue, 2));
pairs = phasePortraitPairs(nStates);
pairScores = nan(size(pairs, 1), 1);
for p = 1:size(pairs, 1)
    predCloud = phasePortraitCloud(xPred(:, 1:nStates), pairs(p, :), options);
    trueCloud = phasePortraitCloud(xTrue(:, 1:nStates), pairs(p, :), options);
    if size(predCloud, 1) >= 2 && size(trueCloud, 1) >= 2
        pairScores(p) = slicedWasserstein2D(predCloud, trueCloud, options.NumProjections, options.TrimFraction);
    else
        pairScores(p) = Inf;
    end
end
finiteScores = pairScores(isfinite(pairScores));
if isempty(finiteScores)
    distanceValue = Inf;
else
    distanceValue = mean(finiteScores);
end
end

function pairs = phasePortraitPairs(nStates)
if nStates >= 2
    pairs = nchoosek(1:nStates, 2);
else
    pairs = [1 1];
end
end

function cloud = phasePortraitCloud(x, pair, options)
nRows = size(x, 1);
if nRows < 2
    cloud = zeros(0, 2);
    return
end
firstRow = 1 + floor(max(0, min(0.9, options.TransientFraction)) * nRows);
firstRow = min(max(firstRow, 1), nRows - 1);
x = x(firstRow:end, :);
x = x(1:max(1, round(options.Subsample)):end, :);
if pair(1) == pair(2)
    cloud = [x(1:end-1, pair(1)), x(2:end, pair(1))];
else
    cloud = x(:, pair);
end
cloud = cloud(all(isfinite(cloud), 2), :);
if size(cloud, 1) > options.MaxPoints
    idx = unique(round(linspace(1, size(cloud, 1), options.MaxPoints)));
    cloud = cloud(idx, :);
end
end

function d = slicedWasserstein2D(P, Q, numProjections, trimFraction)
if size(P, 1) < 2 || size(Q, 1) < 2
    d = Inf;
    return
end
theta = ((0:(numProjections-1)).' + 0.5) * pi / numProjections;
directions = [cos(theta), sin(theta)];
vals = nan(numProjections, 1);
for k = 1:numProjections
    p = sort(P * directions(k, :).');
    q = sort(Q * directions(k, :).');
    if numel(p) ~= numel(q)
        nInterp = max(numel(p), numel(q));
        t = linspace(0, 1, nInterp).';
        p = interp1(linspace(0, 1, numel(p)).', p, t, 'linear', 'extrap');
        q = interp1(linspace(0, 1, numel(q)).', q, t, 'linear', 'extrap');
    end
    vals(k) = mean((p - q).^2, 'omitnan');
end
vals = sort(vals(isfinite(vals)));
if isempty(vals)
    d = Inf;
    return
end
nTrim = floor(max(0, min(0.45, trimFraction)) * numel(vals));
if 2*nTrim < numel(vals)
    vals = vals((1+nTrim):(end-nTrim));
end
d = sqrt(mean(vals, 'omitnan'));
end

function plotSystemLossFigure(taskName, lossHistories, figureDir)
fig = figure('Color', 'w');
ax = axes(fig);
hold(ax, 'on');
colors = lines(max(1, numel(lossHistories)));
hasTraining = false;
hasValidation = false;
legendHandles = gobjects(0);
legendLabels = {};
for k = 1:numel(lossHistories)
    H = lossHistories{k};
    if isempty(H)
        continue
    end
    color = colors(k, :);
    if ~isempty(H.trainingEpoch) && ~isempty(H.trainingLoss)
        hTrain = plot(ax, H.trainingEpoch, positiveForLogPlot(H.trainingLoss), '-', ...
            'Color', color, 'LineWidth', 0.8);
        hasTraining = true;
        legendHandles(end+1) = hTrain; %#ok<AGROW>
        legendLabels{end+1} = sprintf('seed %d train', H.seed); %#ok<AGROW>
    end
    if ~isempty(H.validationEpoch) && ~isempty(H.validationLoss)
        hVal = plot(ax, H.validationEpoch, positiveForLogPlot(H.validationLoss), 'o', ...
            'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 3, 'LineStyle', 'none');
        hasValidation = true;
        legendHandles(end+1) = hVal; %#ok<AGROW>
        legendLabels{end+1} = sprintf('seed %d val', H.seed); %#ok<AGROW>
    end
end
hold(ax, 'off');
set(ax, 'YScale', 'log');
grid(ax, 'on');
xlabel(ax, 'Epoch');
ylabel(ax, 'Loss');
title(ax, sprintf('%s derivative-field training and validation loss', taskName));
if ~isempty(legendHandles)
    legend(ax, legendHandles, legendLabels, 'Location', 'best');
end
if ~hasTraining && ~hasValidation
    text(ax, 0.5, 0.5, 'No saved loss history found', 'HorizontalAlignment', 'center', 'Units', 'normalized');
end
safeTask = matlab.lang.makeValidName(char(taskName));
saveas(fig, fullfile(figureDir, sprintf('%s_training_validation_loss.png', safeTask)));
end

function y = positiveForLogPlot(y)
y = double(y);
finitePositive = y(isfinite(y) & y > 0);
if isempty(finitePositive)
    replacement = eps;
else
    replacement = min(finitePositive) * 0.1;
end
y(~isfinite(y) | y <= 0) = replacement;
end

function plotTrajectoryFigure(t, xTrue, xPred, dynamicsName, figureDir, safeName, wd)
fig = figure('Color', 'w');
tiledlayout(fig, size(xTrue, 2), 1, 'TileSpacing', 'compact');
for j = 1:size(xTrue, 2)
    nexttile;
    plot(t, xTrue(:, j), 'k', 'LineWidth', 1); hold on;
    plot(t, xPred(:, j), 'r--', 'LineWidth', 1); hold off;
    ylabel(sprintf('x_%d', j));
    if j == 1
        title(sprintf('%s derivative-field closed-loop trajectory, WD %.4g', dynamicsName, wd));
    end
    if j == size(xTrue, 2), xlabel('Time'); end
    grid on;
    legend('True', 'Closed-loop', 'Location', 'best');
end
saveas(fig, fullfile(figureDir, sprintf('%s_closed_loop_trajectory.png', safeName)));
end

function plotPhaseFigure(xTrue, xPred, dynamicsName, figureDir, safeName, wd)
nStates = size(xTrue, 2);
fig = figure('Color', 'w');
if nStates >= 2
    pairs = nchoosek(1:nStates, 2);
    tiledlayout(fig, 1, size(pairs, 1), 'TileSpacing', 'compact');
    for p = 1:size(pairs, 1)
        nexttile;
        a = pairs(p, 1); b = pairs(p, 2);
        plot(xTrue(:, a), xTrue(:, b), 'k', 'LineWidth', 1); hold on;
        plot(xPred(:, a), xPred(:, b), 'r--', 'LineWidth', 1); hold off;
        xlabel(sprintf('x_%d', a)); ylabel(sprintf('x_%d', b));
        title(sprintf('x_%d vs x_%d', a, b));
        grid on; axis tight;
        legend('True', 'Closed-loop', 'Location', 'best');
    end
else
    plot(xTrue(:, 1), xPred(:, 1), 'k.', 'MarkerSize', 4);
    xlabel('True x_1'); ylabel('Closed-loop x_1');
    title('x_1 prediction plane');
    grid on; axis tight;
end
sgtitle(sprintf('%s 2D closed-loop phase portraits, WD %.4g', dynamicsName, wd));
saveas(fig, fullfile(figureDir, sprintf('%s_closed_loop_phase_2d.png', safeName)));
end

function plotReturnMapFigure(xTrue, xPred, dynamicsName, figureDir, safeName, wd)
nStates = size(xTrue, 2);
fig = figure('Color', 'w');
tiledlayout(fig, 1, nStates, 'TileSpacing', 'compact');
for j = 1:nStates
    nexttile;
    truePeaks = localPeakValues(xTrue(:, j));
    predPeaks = localPeakValues(xPred(:, j));
    h = gobjects(0);
    labels = {};
    if numel(truePeaks) >= 2
        h(end+1) = plot(truePeaks(1:end-1), truePeaks(2:end), 'k.', 'MarkerSize', 8); hold on; %#ok<AGROW>
        labels{end+1} = 'True'; %#ok<AGROW>
    end
    if numel(predPeaks) >= 2
        h(end+1) = plot(predPeaks(1:end-1), predPeaks(2:end), 'r.', 'MarkerSize', 8); hold on; %#ok<AGROW>
        labels{end+1} = 'Closed-loop'; %#ok<AGROW>
    end
    if ~isempty(h)
        hold off;
        legend(h, labels, 'Location', 'best');
    else
        text(0.5, 0.5, 'Fewer than two peaks', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
    xlabel(sprintf('x_%d peak n', j)); ylabel(sprintf('x_%d peak n+1', j));
    title(sprintf('x_%d peak return', j));
    grid on; axis tight;
end
sgtitle(sprintf('%s peak-to-peak 2D closed-loop return maps, WD %.4g', dynamicsName, wd));
saveas(fig, fullfile(figureDir, sprintf('%s_peak_return_maps_2d.png', safeName)));
end

function peaks = localPeakValues(x)
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    peaks = [];
    return
end
peakMask = x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end);
peaks = x(find(peakMask) + 1);
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
