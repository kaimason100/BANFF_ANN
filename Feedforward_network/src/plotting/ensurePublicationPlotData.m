function plotDataDir = ensurePublicationPlotData(repoRoot, plotDataDir, options)
% ensurePublicationPlotData  Build legacy publication-plot MAT files.
%
% The publication figure scripts pre-date the package layout and expect files
% such as MNIST_plot_data.mat in their own folder. This helper recreates those
% files from the current saved-network/data layout without changing the figure
% scripts' plotting commands.

arguments
    repoRoot (1,:) char
    plotDataDir (1,:) char
    options.UseSeededNetworks (1,1) logical = true
    options.Seed (1,1) double = 0
    options.Force (1,1) logical = false
    options.ClosedLoopTspan (1,2) double = [0 5000]
    options.LorenzClosedLoopTspan (1,2) double = [0 10000]
    options.ClosedLoopOutputDt = []
    options.ManualInitialCondition = []
    options.DynamicalSystemNetworkSet (1,1) string = "seeded_state_random_derivative_ode45_lr_0p01"
    options.UseContinuationIfAvailable (1,1) logical = true
    options.ContinuationNetworkSets = "auto"
    options.TestICRandomSeed (1,1) double = 9100
    options.TestICPerturbationScale (1,1) double = 0
    options.MnistBatchSize (1,1) double = 256
    options.Tasks = "all"
end

if ~exist(plotDataDir, 'dir'), mkdir(plotDataDir); end
addpath(fullfile(repoRoot, 'src', 'preprocessing'));
addpath(fullfile(repoRoot, 'src', 'initializers'));
addpath(fullfile(repoRoot, 'src', 'plotting'));
addpath(fullfile(repoRoot, 'src', 'dynamical_systems'));

tasks = normalizeTasks(options.Tasks);
runAll = any(strcmp(tasks, "all"));
runMaybe(runAll || any(strcmp(tasks, "mnist")), @() generateMnistPlotData(repoRoot, plotDataDir, options), 'MNIST', runAll);
runMaybe(runAll || any(strcmp(tasks, "toyota")), @() generateToyotaPlotData(repoRoot, plotDataDir, options), 'Toyota', runAll);
runMaybe(runAll || any(strcmp(tasks, "dynamical_systems")), @() generateDynamicalSystemPlotData(repoRoot, plotDataDir, options), 'dynamical systems', runAll);
runMaybe(runAll || any(strcmp(tasks, "pong")), @() generatePongPlotData(repoRoot, plotDataDir, options), 'Pong', runAll);
runMaybe(runAll || any(strcmp(tasks, "lqr")), @() generateLqrPlotData(repoRoot, plotDataDir, options), 'LQR', runAll);
runMaybe(runAll || any(strcmp(tasks, "bias_histogram")), @() generateBiasHistogramPlotData(repoRoot, plotDataDir, options), 'bias histogram', runAll);
end

function tasks = normalizeTasks(tasks)
if ischar(tasks), tasks = string({tasks}); end
tasks = lower(string(tasks));
tasks = strrep(tasks, "-", "_");
tasks(tasks == "motor_control") = "lqr";
tasks(tasks == "ds") = "dynamical_systems";
tasks(tasks == "dynamics") = "dynamical_systems";
tasks(tasks == "bias") = "bias_histogram";
end

function runMaybe(shouldRun, generatorFcn, label, continueOnError)
if ~shouldRun, return; end
try
    generatorFcn();
catch ME
    if continueOnError
        warningId = ['PublicationPlotData:' matlab.lang.makeValidName(label)];
        warning(warningId, 'Could not generate %s publication plot data: %s', label, ME.message);
    else
        rethrow(ME);
    end
end
end

function generateMnistPlotData(repoRoot, plotDataDir, options)
outFile = fullfile(plotDataDir, 'MNIST_plot_data.mat');
if isfile(outFile) && ~options.Force, return; end
modelPath = resolveModelPath(repoRoot, 'mnist', 'mnist_network.mat', options);
Snet = load(modelPath);
net = loadNetworkObject(Snet, modelPath);
D = load(fullfile(repoRoot, 'data', 'mnist.mat'));
[~, ~, XTest, YTest, numOut] = mnistData(D.training, D.test);
if isfield(Snet, 'imageStats')
    XEval = applyImageStandardization(XTest, Snet.imageStats);
else
    XEval = XTest;
end
scores = numericOutput(predictImageBatches(net, XEval, options.MnistBatchSize));
plotData.output = scoresToClassIndex(scores, numOut, numel(YTest));
plotData.true = double(YTest(:));
exampleRows = firstClassExamples(plotData.true, numOut);
plotData.activationExampleRows = exampleRows;
plotData.activationExamples = hiddenActivations(net, subsetImages(XEval, exampleRows), "last");
save(outFile, 'plotData', '-v7.3');
end

function generateToyotaPlotData(repoRoot, plotDataDir, options)
outFile = fullfile(plotDataDir, 'Toyota_plot_data.mat');
if isfile(outFile) && ~options.Force, return; end
modelPath = resolveModelPath(repoRoot, 'toyota', 'Toyota_network.mat', options);
Snet = load(modelPath);
net = loadNetworkObject(Snet, modelPath);
[X, Y] = loadRegressionData(repoRoot, 'toyota');
assert(isfield(Snet, 'split') && isfield(Snet.split, 'testIdx'), 'Toyota network lacks test split metadata.');
XTest = applyFeatureStandardization(X(:, Snet.split.testIdx), Snet.featureStats);
YPredNorm = numericOutput(predictFeatureBatch(net, XTest));
plotData.output = inverseTargetStandardization(YPredNorm(:), Snet.targetStats);
plotData.true = Y(Snet.split.testIdx).';
plotData.activations = hiddenActivations(net, XTest, "last");
save(outFile, 'plotData', '-v7.3');
end

function generateDynamicalSystemPlotData(repoRoot, plotDataDir, options)
names = ["Lorenz" "SprottB" "SprottC" "SprottS" "Rikitake" "MO0" "MO5" "MO7" "MO13"];
for name = names
    outFile = fullfile(plotDataDir, sprintf('%s_plot_data.mat', char(name)));
    videoFile = fullfile(plotDataDir, sprintf('%s_video_plot_data.mat', char(name)));
    needsVideoFile = name == "Lorenz";
    if isfile(outFile) && (~needsVideoFile || isfile(videoFile)) && ~options.Force, continue; end
    dynamicsOptions = options;
    if name == "Lorenz"
        dynamicsOptions.ClosedLoopTspan = options.LorenzClosedLoopTspan;
    end
    [modelPath, selectedNetworkSet, usedContinuation] = resolveDynamicalSystemModelPath(repoRoot, char(name), dynamicsOptions);
    Snet = load(modelPath);
    net = loadNetworkObject(Snet, modelPath);
    [t, xTrue, xPred, activations1, activations2] = dynamicalClosedLoop(repoRoot, net, Snet, char(name), dynamicsOptions);
    plotData.time = t;
    plotData.true = xTrue;
    plotData.output = xPred;
    plotData.activations = activations2;
    plotData.activations_hidden1 = activations1;
    plotData.activations_hidden2 = activations2;
    [plotData.biases_hidden1, plotData.biases_hidden2] = hiddenBiases(net);
    plotData.method = 'state-random-derivative-ode45';
    plotData.networkSet = char(selectedNetworkSet);
    plotData.usedContinuation = usedContinuation;
    plotData.networkPath = modelPath;
    plotData.closedLoopTspan = dynamicsOptions.ClosedLoopTspan;
    plotData.closedLoopOutputDt = dynamicsOptions.ClosedLoopOutputDt;
    plotData.testICRandomSeed = dynamicsOptions.TestICRandomSeed;
    plotData.testICPerturbationScale = dynamicsOptions.TestICPerturbationScale;
    save(outFile, 'plotData', '-v7.3');
    if name == "Lorenz"
        save(videoFile, 'plotData', '-v7.3');
    end
end
end

function generatePongPlotData(repoRoot, plotDataDir, options)
outFile = fullfile(plotDataDir, 'Pong_plot_data.mat');
videoFile = fullfile(plotDataDir, 'Pong_video_plot_data.mat');
if isfile(outFile) && isfile(videoFile) && ~options.Force, return; end
modelPath = resolveModelPath(repoRoot, 'pong', 'Pong_network.mat', options);
Snet = load(modelPath);
net = loadNetworkObject(Snet, modelPath);
[plotData, videoPlotData] = simulatePongForPlot(net);
plotData.biases_hidden1 = videoPlotData.biases_hidden1;
plotData.biases_hidden2 = videoPlotData.biases_hidden2;
save(outFile, 'plotData', '-v7.3');
plotData = videoPlotData; %#ok<NASGU>
save(videoFile, 'plotData', '-v7.3');
end

function generateLqrPlotData(repoRoot, plotDataDir, options)
outFile = fullfile(plotDataDir, 'LQR_plot_data.mat');
videoFile = fullfile(plotDataDir, 'LQR_video_plot_data.mat');
if isfile(outFile) && isfile(videoFile) && ~options.Force, return; end
modelPath = resolveModelPath(repoRoot, 'motor_control', 'LQR_network.mat', options);
Snet = load(modelPath);
net = loadNetworkObject(Snet, modelPath);
target = selectLqrPlotTarget(Snet.testSequence);
sim = simulateMotorControllerForPlot(net, Snet, target);
plotData.true = sim.eeLqr.';
plotData.output = sim.eeNet.';
plotData.activations = sim.activations;
plotData.activations_hidden1 = sim.activations_hidden1;
plotData.activations_hidden2 = sim.activations_hidden2;
plotData.linkLengths = [sim.l1 sim.l2];
plotData.t = sim.t;
plotData.dt = sim.dt;
plotData.lqrTrajectory = sim.xLqr;
plotData.netTrajectory = sim.xNet;
plotData.targetEE = {target};
[plotData.biases_hidden1, plotData.biases_hidden2] = hiddenBiases(net);
save(outFile, 'plotData', '-v7.3');
save(videoFile, 'plotData', '-v7.3');
end

function target = selectLqrPlotTarget(testSequence)
if isfield(testSequence, 'targetEE')
    if iscell(testSequence.targetEE)
        target = testSequence.targetEE{min(2, numel(testSequence.targetEE))};
    else
        target = testSequence.targetEE(min(2, size(testSequence.targetEE, 1)), :).';
    end
elseif isfield(testSequence, 'testTargets')
    target = testSequence.testTargets(min(2, size(testSequence.testTargets, 1)), :).';
else
    error('LQR testSequence lacks targetEE/testTargets for publication plot data.');
end
target = target(:);
end

function pathOut = resolveModelPath(repoRoot, task, singleName, options)
if options.UseSeededNetworks
    taskNames = seededTaskAliases(task);
    candidates = {};
    for i = 1:numel(taskNames)
        taskName = taskNames{i};
        candidates{end+1} = fullfile(repoRoot, 'trained_networks', 'seeded_closed_loop_random_windows', taskName, sprintf('%s_seed_%03d_network.mat', taskName, options.Seed)); %#ok<AGROW>
        candidates{end+1} = fullfile(repoRoot, 'trained_networks', 'seeded', taskName, sprintf('%s_seed_%03d_network.mat', taskName, options.Seed)); %#ok<AGROW>
        candidates{end+1} = fullfile(repoRoot, 'trained_networks', 'seeded_closed_loop', taskName, sprintf('%s_seed_%03d_network.mat', taskName, options.Seed)); %#ok<AGROW>
    end
else
    candidates = {fullfile(repoRoot, 'trained_networks', singleName)};
end
existsMask = cellfun(@isfile, candidates);
idx = find(existsMask, 1, 'first');
assert(~isempty(idx), 'Missing saved network for %s. Checked: %s', task, strjoin(candidates, ', '));
pathOut = candidates{idx};
end

function [modelPath, selectedNetworkSet, usedContinuation] = resolveDynamicalSystemModelPath(repoRoot, task, options)
if ~options.UseSeededNetworks
    selectedNetworkSet = "";
    usedContinuation = false;
    modelPath = fullfile(repoRoot, 'trained_networks', sprintf('%s_network.mat', char(task)));
    assert(isfile(modelPath), 'Missing saved network for %s: %s', task, modelPath);
    return
end
seedValue = double(options.Seed);
baseNetworkSet = string(options.DynamicalSystemNetworkSet);
usedContinuation = false;
selectedNetworkSet = baseNetworkSet;
if options.UseContinuationIfAvailable
    continuationSets = resolveDynamicalContinuationSets(repoRoot, options.ContinuationNetworkSets);
    for k = 1:numel(continuationSets)
        candidate = seededDynamicalSystemPath(repoRoot, continuationSets(k), task, seedValue);
        if isfile(candidate)
            modelPath = candidate;
            selectedNetworkSet = continuationSets(k);
            usedContinuation = true;
            fprintf('%s publication plot data: using continuation network set %s\n', task, char(selectedNetworkSet));
            return
        end
    end
end
modelPath = seededDynamicalSystemPath(repoRoot, baseNetworkSet, task, seedValue);
assert(isfile(modelPath), 'Missing derivative-field saved network for %s seed %d: %s', task, seedValue, modelPath);
fprintf('%s publication plot data: using network set %s\n', task, char(selectedNetworkSet));
end

function modelPath = seededDynamicalSystemPath(repoRoot, networkSet, task, seedValue)
folderName = matlab.lang.makeValidName(char(task));
modelPath = fullfile(repoRoot, 'trained_networks', char(networkSet), folderName, ...
    sprintf('%s_seed_%03d_network.mat', folderName, double(seedValue)));
end

function candidateSets = resolveDynamicalContinuationSets(repoRoot, continuationNetworkSets)
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

function names = seededTaskAliases(task)
task = char(task);
switch lower(task)
    case 'pong'
        names = {'Pong', 'pong'};
    case {'motor_control', 'lqr', 'lqr_two_link_arm'}
        names = {'LQR_two_link_arm', 'motor_control', 'LQR'};
    otherwise
        valid = matlab.lang.makeValidName(task);
        names = unique({task, lower(task), valid}, 'stable');
end
end

function net = loadNetworkObject(S, modelPath)
if isfield(S, 'net')
    net = S.net;
elseif isfield(S, 'dlnet')
    net = S.dlnet;
elseif isfield(S, 'trainedNet')
    net = S.trainedNet;
elseif isfield(S, 'network')
    net = S.network;
elseif isfield(S, 'netObj')
    net = S.netObj;
elseif isfield(S, 'bestNet')
    net = S.bestNet;
else
    error('%s has no loadable network object.', modelPath);
end
end

function [XTrain, YTrain, XTest, YTest, numOut] = mnistData(training, test)
XTrain = ensure4d(training.images);
XTest = ensure4d(test.images);
YTrain = categorical(double(training.labels(:)) + 1);
YTest = categorical(double(test.labels(:)) + 1);
numOut = numel(categories(YTrain));
end

function X = ensure4d(X)
if ndims(X) == 3
    X = reshape(X, size(X,1), size(X,2), 1, size(X,3));
end
X = single(X);
end

function XSub = subsetImages(X, rows)
XSub = X(:, :, :, rows);
end

function rows = firstClassExamples(y, numOut)
rows = zeros(numOut, 1);
for i = 1:numOut
    idx = find(y == i, 1, 'first');
    assert(~isempty(idx), 'MNIST class %d has no example.', i);
    rows(i) = idx;
end
end

function XOut = applyImageStandardization(X, stats)
scale = 1;
if isfield(stats, 'scale'), scale = stats.scale; end
XOut = scale * ((single(X) - single(stats.mu)) ./ single(stats.sigma));
end

function XOut = applyFeatureStandardization(X, stats)
scale = 1;
if isfield(stats, 'scale'), scale = stats.scale; end
XOut = scale * ((X - stats.mu) ./ stats.sigma);
end

function Y = inverseTargetStandardization(YNorm, stats)
if isfield(stats, 'sigma')
    sigma = stats.sigma;
elseif isfield(stats, 'yStd')
    sigma = stats.yStd;
else
    error('Target standardization metadata lacks sigma/yStd.');
end
if isfield(stats, 'mu')
    mu = stats.mu;
elseif isfield(stats, 'yMu')
    mu = stats.yMu;
else
    error('Target standardization metadata lacks mu/yMu.');
end
Y = YNorm .* sigma + mu;
end

function [X, Y] = loadRegressionData(repoRoot, task)
switch lower(task)
    case 'toyota'
        S = load(fullfile(repoRoot, 'data', 'toyota_dataset.mat'));
        D = S.data;
        D(~isfinite(D)) = 0;
        Y = D(:, 3);
        X = D(:, [1, 2, 4:size(D, 2)]).';
    otherwise
        error('Unsupported regression task for publication plots: %s', task);
end
end

function Y = predictImageBatches(net, X, batchSize)
n = size(X, 4);
Y = [];
for i = 1:batchSize:n
    idx = i:min(n, i+batchSize-1);
    XBatch = X(:, :, :, idx);
    if networkUsesImageInput(net)
        YBatch = numericOutput(predict(net, single(XBatch)));
    else
        XFlat = reshape(single(XBatch), [], numel(idx));
        YBatch = numericOutput(predictFeatureBatch(net, XFlat));
    end
    YBatch = orientScoreBatch(YBatch, numel(idx));
    Y = [Y; YBatch]; %#ok<AGROW>
end
end

function scores = orientScoreBatch(scores, nObs)
% Return classification scores as observations-by-classes.
if size(scores, 1) == nObs
    return;
end
if size(scores, 2) == nObs
    scores = scores.';
    return;
end
error('Unexpected batch score size [%s] for %d observations.', num2str(size(scores)), nObs);
end

function Y = predictFeatureBatch(net, X)
assert(isnumeric(X) && ismatrix(X), 'Feature predictors must be features-by-observations.');
if isa(net, 'dlnetwork')
    Y = predict(net, dlarray(single(X), 'CB'));
else
    Y = predict(net, single(X.'));
end
Y = numericOutput(Y);
end

function Y = numericOutput(Y)
if isa(Y, 'dlarray'), Y = extractdata(Y); end
if isa(Y, 'gpuArray'), Y = gather(Y); end
if istable(Y), Y = table2array(Y); end
Y = double(Y);
end

function tf = networkUsesImageInput(net)
tf = false;
if isprop(net, 'Layers') && ~isempty(net.Layers)
    tf = contains(lower(class(net.Layers(1))), 'imageinput');
end
end

function predIdx = scoresToClassIndex(scores, numOut, nObs)
if size(scores, 1) == nObs && size(scores, 2) == numOut
    [~, predIdx] = max(scores, [], 2);
elseif size(scores, 2) == nObs && size(scores, 1) == numOut
    [~, predIdx] = max(scores, [], 1);
    predIdx = predIdx(:);
else
    error('Unexpected score size [%s].', num2str(size(scores)));
end
end

function A = hiddenActivations(net, X, whichLayer)
if isa(net, 'dlnetwork')
    A = dlnetworkHiddenActivations(net, X, whichLayer);
    return
end
layerName = hiddenLayerName(net, whichLayer);
try
    if ndims(X) == 4
        A = activations(net, single(X), layerName, 'OutputAs', 'rows');
    else
        A = activations(net, single(X.'), layerName, 'OutputAs', 'rows');
    end
catch
    A = fallbackHiddenActivations(net, X, whichLayer);
end
A = numericOutput(A);
if size(A, 1) > 0 && size(A, 2) < size(A, 1) && size(A, 2) <= 3
    A = A.';
end
end

function A = dlnetworkHiddenActivations(net, X, whichLayer)
if ndims(X) == 4
    X = reshape(single(X), [], size(X, 4));
else
    X = single(X);
end
learnables = net.Learnables;
weightRows = find(strcmp(string(learnables.Parameter), "Weights"));
biasRows = find(strcmp(string(learnables.Parameter), "Bias"));
assert(numel(weightRows) >= 2 && numel(biasRows) >= 2, 'dlnetwork hidden activations require at least two dense layers.');
W1 = numericOutput(learnables.Value{weightRows(1)});
b1 = numericOutput(learnables.Value{biasRows(1)});
W2 = numericOutput(learnables.Value{weightRows(2)});
b2 = numericOutput(learnables.Value{biasRows(2)});
if size(X, 1) ~= size(W1, 2) && size(X, 2) == size(W1, 2)
    X = X.';
end
assert(size(X, 1) == size(W1, 2), 'dlnetwork activation input has %d features; expected %d.', size(X, 1), size(W1, 2));
A1 = tanh(single(W1) * single(X) + single(b1));
if whichLayer == "first"
    A = A1.';
else
    A = tanh(single(W2) * A1 + single(b2)).';
end
A = numericOutput(A);
end

function layerName = hiddenLayerName(net, whichLayer)
assert(isprop(net, 'Layers'), 'Hidden activations require a SeriesNetwork/DAGNetwork with Layers.');
classes = arrayfun(@(L) lower(class(L)), net.Layers, 'UniformOutput', false);
idx = find(contains(classes, 'tanh'));
if isempty(idx)
    idx = find(contains(classes, 'relu'));
end
assert(~isempty(idx), 'No hidden activation layer found.');
if whichLayer == "first"
    layerName = net.Layers(idx(1)).Name;
else
    layerName = net.Layers(idx(end)).Name;
end
end

function A = fallbackHiddenActivations(net, X, whichLayer)
assert(isprop(net, 'Layers'), 'Fallback activations require Layers.');
if ndims(X) == 4
    X = reshape(single(X), [], size(X, 4));
end
expectedInputs = prod(double(net.Layers(1).InputSize));
if size(X, 1) ~= expectedInputs && size(X, 2) == expectedInputs
    X = X.';
end
fc = find(arrayfun(@(L) contains(lower(class(L)), 'fullyconnected'), net.Layers));
if whichLayer == "first"
    k = fc(1);
else
    k = fc(2);
end
W = net.Layers(k).Weights;
b = net.Layers(k).Bias;
if k == fc(1)
    A = tanh(W * single(X) + b);
else
    A1 = tanh(net.Layers(fc(1)).Weights * single(X) + net.Layers(fc(1)).Bias);
    A = tanh(W * A1 + b);
end
A = A.';
end

function [b1, b2] = hiddenBiases(net)
if isa(net, 'dlnetwork') && isprop(net, 'Learnables')
    learnables = net.Learnables;
    biasRows = strcmp(string(learnables.Parameter), "Bias");
    assert(nnz(biasRows) >= 2, 'dlnetwork has fewer than two Bias learnables.');
    biasValues = learnables.Value(biasRows);
    b1 = numericOutput(biasValues{1});
    b2 = numericOutput(biasValues{2});
else
    fc = find(arrayfun(@(L) contains(lower(class(L)), 'fullyconnected'), net.Layers));
    assert(numel(fc) >= 2, 'Network has fewer than two fully connected layers.');
    b1 = numericOutput(net.Layers(fc(1)).Bias);
    b2 = numericOutput(net.Layers(fc(2)).Bias);
end
b1 = b1(:);
b2 = b2(:);
end

function generateBiasHistogramPlotData(repoRoot, plotDataDir, options)
outFile = fullfile(plotDataDir, 'Bias_histogram_plot_data.mat');
if isfile(outFile) && ~options.Force, return; end
tasks = {
    'mnist', 'mnist_network.mat'
    'toyota', 'Toyota_network.mat'
    'motor_control', 'LQR_network.mat'
    'pong', 'Pong_network.mat'
    'Lorenz', 'Lorenz_network.mat'
};
b = nan(16000, size(tasks, 1));
for i = 1:size(tasks, 1)
    if strcmpi(tasks{i, 1}, 'Lorenz')
        modelPath = resolveDynamicalSystemModelPath(repoRoot, tasks{i, 1}, options);
    else
        modelPath = resolveModelPath(repoRoot, tasks{i, 1}, tasks{i, 2}, options);
    end
    Snet = load(modelPath);
    net = loadNetworkObject(Snet, modelPath);
    [~, b2] = hiddenBiases(net);
    n = min(size(b, 1), numel(b2));
    b(1:n, i) = b2(1:n);
end
save(outFile, 'b', '-v7.3');
end

function [t, xTrue, xPred, A1, A2] = dynamicalClosedLoop(repoRoot, net, Snet, dynamicsName, options)
assert(isfield(Snet, 'stateStats'), '%s saved network lacks stateStats.', dynamicsName);
assert(isa(net, 'dlnetwork'), '%s publication plot data requires the derivative-field dlnetwork.', dynamicsName);
dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
row = find(strcmp(string(dyn{:, 1}), string(dynamicsName)), 1, 'first');
assert(~isempty(row), 'Could not find %s in dynamics list.', dynamicsName);
d = dyn{row, 2};
ic = dyn{row, 4:6}.';
icBase = selectDynamicalInitialCondition(ic, d, options.ManualInitialCondition, dynamicsName);
x0Norm = applyStateStandardization(icBase, Snet.stateStats);
if options.TestICPerturbationScale > 0
    stream = RandStream('mt19937ar', 'Seed', double(options.TestICRandomSeed));
    x0Norm = x0Norm + double(options.TestICPerturbationScale) * (2*rand(stream, size(x0Norm)) - 1);
end
x0Norm = x0Norm(:);

odeOptions = odeset();
tRequest = odeRequestTimes(options.ClosedLoopTspan, options.ClosedLoopOutputDt);
if numel(tRequest) == 2
    [t, xTrue] = ode45(@(tt, xx) trueStandardizedDerivative(dynamicsName, xx, Snet.stateStats), tRequest, x0Norm, odeOptions);
    [tPred, xPred] = ode45(@(tt, xx) learnedDerivativeFieldPredict(net, xx), t, x0Norm, odeOptions);
    assert(isequal(size(tPred), size(t)), '%s learned trajectory did not return on the true trajectory time base.', dynamicsName);
else
    [t, xTrue] = ode45(@(tt, xx) trueStandardizedDerivative(dynamicsName, xx, Snet.stateStats), tRequest, x0Norm, odeOptions);
    [tPred, xPred] = ode45(@(tt, xx) learnedDerivativeFieldPredict(net, xx), tRequest, x0Norm, odeOptions);
    assert(isequal(size(tPred), size(t)), '%s learned and true trajectories have different output times.', dynamicsName);
end
assert(isequal(size(xPred), size(xTrue)), '%s learned and true trajectories differ in size.', dynamicsName);
A1 = hiddenActivations(net, xPred.', "first");
A2 = hiddenActivations(net, xPred.', "last");
end

function icBase = selectDynamicalInitialCondition(defaultIc, stateDim, manualInitialCondition, dynamicsName)
if isempty(manualInitialCondition)
    icBase = defaultIc(1:stateDim).';
    return;
end
assert(isnumeric(manualInitialCondition) && isvector(manualInitialCondition), '%s ManualInitialCondition must be a numeric vector.', dynamicsName);
assert(numel(manualInitialCondition) >= stateDim, '%s ManualInitialCondition must contain at least %d entries.', dynamicsName, stateDim);
icBase = manualInitialCondition(1:stateDim);
icBase = icBase(:).';
assert(all(isfinite(icBase)), '%s ManualInitialCondition contains non-finite values.', dynamicsName);
end

function xOut = applyStateStandardization(x, stats)
scale = 1;
if isfield(stats, 'scale'), scale = stats.scale; end
xOut = scale * ((double(x) - stats.mu) ./ stats.sigma);
xOut(~isfinite(xOut)) = 0;
end

function xRaw = inverseStateStandardization(xNorm, stats)
scale = 1;
if isfield(stats, 'scale'), scale = stats.scale; end
xRaw = (double(xNorm) ./ scale) .* stats.sigma + stats.mu;
end

function tRequest = odeRequestTimes(tspan, outputDt)
assert(isnumeric(tspan) && numel(tspan) == 2 && tspan(2) > tspan(1), 'ClosedLoopTspan must be [t0 tf].');
if isempty(outputDt)
    tRequest = double(tspan(:));
    return
end
assert(isscalar(outputDt) && isfinite(outputDt) && outputDt > 0, 'ClosedLoopOutputDt must be empty or a positive scalar.');
nWholeSteps = floor((double(tspan(2)) - double(tspan(1))) / double(outputDt));
tRequest = double(tspan(1)) + (0:nWholeSteps).' * double(outputDt);
if tRequest(end) < double(tspan(2))
    tRequest(end + 1, 1) = double(tspan(2));
end
end

function dxNormDt = trueStandardizedDerivative(dynamicsName, xNorm, stats)
xNorm = double(xNorm(:));
xRaw = inverseStateStandardization(xNorm.', stats).';
dxRawDt = dynamicsDerivative(dynamicsName, xRaw);
scale = 1;
if isfield(stats, 'scale'), scale = stats.scale; end
dxNormDt = scale * (dxRawDt(:) ./ stats.sigma(:));
assert(all(isfinite(dxNormDt)), '%s true standardized derivative produced non-finite values.', dynamicsName);
end

function dxdt = learnedDerivativeFieldPredict(net, x)
x = double(x(:));
y = numericOutput(predict(net, dlarray(single(x), 'CB')));
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

function [plotData, videoPlotData] = simulatePongForPlot(net)
rng(0, 'twister');
dt = 0.01;
nSteps = 3600;
ballSpeed = 0.02;
maxBounceAngle = pi/4;
paddleHeight = 0.2;
paddleWidth = 0.02;
ballRadius = 0.02;
oppPaddleX = 0.05;
netPaddleX = 1 - oppPaddleX - paddleWidth;
ball = [0.5 0.5];
theta = (pi/4)*rand() + pi/8;
vel = ballSpeed * [cos(theta), sin(theta)];
oppY = 0.5;
netY = 0.5;
paddleSpeed = 0.025;
ballSeq = zeros(nSteps, 2);
oppSeq = zeros(nSteps, 1);
netSeq = zeros(nSteps, 1);
states = zeros(nSteps, 5);
actions = strings(nSteps, 1);
networkHits = 0;
networkMisses = 0;
opponentMisses = 0;
for k = 1:nSteps
    state = [ball vel netY];
    states(k,:) = state;
    action = predictPongAction(net, state);
    actions(k) = action;
    if action == "Up", netY = max(paddleHeight/2, netY - paddleSpeed); end
    if action == "Down", netY = min(1-paddleHeight/2, netY + paddleSpeed); end
    oppY = min(max(ball(2), paddleHeight/2), 1-paddleHeight/2);
    ball = ball + vel;
    if ball(2) - ballRadius <= 0
        ball(2) = ballRadius;
        vel(2) = abs(vel(2));
    elseif ball(2) + ballRadius >= 1
        ball(2) = 1 - ballRadius;
        vel(2) = -abs(vel(2));
    end
    if vel(1) < 0 && ball(1) - ballRadius <= oppPaddleX + paddleWidth
        if ball(2) >= oppY - paddleHeight/1.5 && ball(2) <= oppY + paddleHeight/1.5
            relIntersect = (ball(2) - oppY) / (paddleHeight/2);
            angle = relIntersect * maxBounceAngle;
            vel(1) = abs(ballSpeed * cos(angle));
            vel(2) = ballSpeed * sin(angle);
            ball(1) = oppPaddleX + paddleWidth + ballRadius;
        else
            opponentMisses = opponentMisses + 1;
            [ball, vel, oppY, netY] = resetPongRound(ballSpeed, paddleHeight);
        end
    end
    if vel(1) > 0 && ball(1) + ballRadius >= netPaddleX
        if ball(2) >= netY - paddleHeight/1.5 && ball(2) <= netY + paddleHeight/1.5
            networkHits = networkHits + 1;
            relIntersect = (ball(2) - netY) / (paddleHeight/2);
            angle = relIntersect * maxBounceAngle;
            vel(1) = -abs(ballSpeed * cos(angle));
            vel(2) = ballSpeed * sin(angle);
            ball(1) = netPaddleX - ballRadius;
        else
            networkMisses = networkMisses + 1;
            [ball, vel, oppY, netY] = resetPongRound(ballSpeed, paddleHeight);
        end
    end
    ball = min(max(ball, [0 0]), [1 1]);
    ballSeq(k,:) = ball;
    oppSeq(k) = oppY;
    netSeq(k) = netY;
end
plotData.true = ballSeq;
plotData.output = netSeq;
plotData.activations = hiddenActivations(net, states.', "last");
videoPlotData = plotData;
videoPlotData.dt = dt;
videoPlotData.paddleHeight = paddleHeight;
videoPlotData.paddleWidth = paddleWidth;
videoPlotData.ballRadius = ballRadius;
videoPlotData.oppPaddleX = oppPaddleX;
videoPlotData.netPaddleX = netPaddleX;
videoPlotData.ballPosSeq = ballSeq;
videoPlotData.oppPaddleYSeq = oppSeq;
videoPlotData.netPaddleYSeq = netSeq;
videoPlotData.actions = actions;
videoPlotData.networkHits = networkHits;
videoPlotData.networkMisses = networkMisses;
videoPlotData.opponentMisses = opponentMisses;
videoPlotData.activations_hidden1 = hiddenActivations(net, states.', "first");
videoPlotData.activations_hidden2 = plotData.activations;
[videoPlotData.biases_hidden1, videoPlotData.biases_hidden2] = hiddenBiases(net);
end

function [ball, vel, oppY, netY] = resetPongRound(ballSpeed, paddleHeight)
ball = [0.5, 0.5];
theta = (pi/4)*rand() + pi/8;
if rand() < 0.5
    vel = ballSpeed * [cos(theta), sin(theta)];
else
    vel = ballSpeed * [-cos(theta), sin(theta)];
end
oppY = 0.5;
netY = 0.5;
oppY = min(max(oppY, paddleHeight/2), 1-paddleHeight/2);
netY = min(max(netY, paddleHeight/2), 1-paddleHeight/2);
end

function action = predictPongAction(net, state)
raw = predictFeatureBatch(net, state(:));
if iscategorical(raw)
    action = string(raw);
    return;
end
[~, idx] = max(raw(:));
labels = ["Up" "Down" "Stay"];
action = labels(idx);
end

function sim = simulateMotorControllerForPlot(net, S, target)
seq = S.testSequence;
nT = motorSequenceLength(seq);
sim.dt = seq.dt;
if isfield(seq, 't')
    sim.t = seq.t(:).';
else
    sim.t = (0:nT-1) * seq.dt;
end
sim.l1 = seq.l1; sim.l2 = seq.l2;
xTargetJoint = inverseKinematicsTarget(target, seq.l1, seq.l2);
xNet = zeros(4, nT); xLqr = zeros(4, nT);
eeNet = zeros(2, nT); eeLqr = zeros(2, nT);
numInputs = motorNumInputs(S, seq);
features = zeros(nT-1, numInputs);
for k = 1:nT-1
    eNet = xNet(:, k) - xTargetJoint;
    eLqr = xLqr(:, k) - xTargetJoint;
    features(k,:) = encodeJointError(eNet, numInputs).';
    uNet = predictMotorControlForPlot(net, features(k,:).', S, numInputs);
    uLqr = -seq.K * eLqr;
    xNet(:, k+1) = xNet(:, k) + seq.dt * (seq.A * xNet(:, k) + seq.B * uNet);
    xLqr(:, k+1) = xLqr(:, k) + seq.dt * (seq.A * xLqr(:, k) + seq.B * uLqr);
    eeNet(:, k) = forwardKinematics(xNet(:, k), seq.l1, seq.l2);
    eeLqr(:, k) = forwardKinematics(xLqr(:, k), seq.l1, seq.l2);
end
eeNet(:, nT) = forwardKinematics(xNet(:, nT), seq.l1, seq.l2);
eeLqr(:, nT) = forwardKinematics(xLqr(:, nT), seq.l1, seq.l2);
sim.xNet = xNet; sim.xLqr = xLqr; sim.eeNet = eeNet; sim.eeLqr = eeLqr;
sim.activations_hidden1 = hiddenActivations(net, features.', "first");
sim.activations_hidden2 = hiddenActivations(net, features.', "last");
sim.activations = sim.activations_hidden2;
end

function nT = motorSequenceLength(seq)
if isfield(seq, 'nT')
    nT = seq.nT;
elseif isfield(seq, 't')
    nT = numel(seq.t);
else
    error('LQR testSequence lacks nT/t for publication plot data.');
end
end

function numInputs = motorNumInputs(S, seq)
if isfield(seq, 'numInputs')
    numInputs = seq.numInputs;
elseif isfield(S, 'mu_X')
    numInputs = numel(S.mu_X);
elseif isfield(seq, 'mu_X')
    numInputs = numel(seq.mu_X);
else
    error('LQR saved network lacks numInputs/mu_X metadata.');
end
end

function u = predictMotorControlForPlot(net, phi, S, numInputs)
if isfield(S, 'mu_X')
    muX = S.mu_X;
    sigmaX = S.sigma_X;
    muY = S.mu_Y;
    sigmaY = S.sigma_Y;
elseif isfield(S, 'testSequence')
    muX = S.testSequence.mu_X;
    sigmaX = S.testSequence.sigma_X;
    muY = S.testSequence.mu_Y;
    sigmaY = S.testSequence.sigma_Y;
else
    error('LQR saved network lacks input/output standardization metadata.');
end
x = (phi - muX) ./ (sigmaX * sqrt(numInputs));
x = reshape(x, 1, numInputs);
if isa(net, 'dlnetwork')
    y = predict(net, dlarray(single(x.'), 'CB'));
else
    y = predict(net, single(x));
end
y = numericOutput(y);
u = y(:) .* sigmaY + muY;
end

function phi = encodeJointError(e, numInputs)
if numInputs == 6
    phi = [sin(e(1)); cos(e(1)); sin(e(2)); cos(e(2)); e(3); e(4)];
elseif numInputs == 8
    phi = [e; sin(e(1)); cos(e(1)); sin(e(2)); cos(e(2))];
else
    error('Unsupported LQR feature count %d for publication plot data.', numInputs);
end
end

function xTarget = inverseKinematicsTarget(target, l1, l2)
x = target(1); y = target(2);
c2 = (x^2 + y^2 - l1^2 - l2^2) / (2*l1*l2);
c2 = min(max(c2, -1), 1);
q2 = acos(c2);
q1 = atan2(y, x) - atan2(l2*sin(q2), l1 + l2*cos(q2));
xTarget = [q1; q2; 0; 0];
end

function ee = forwardKinematics(x, l1, l2)
q1 = x(1); q2 = x(2);
ee = [l1*cos(q1) + l2*cos(q1+q2); l1*sin(q1) + l2*sin(q1+q2)];
end
