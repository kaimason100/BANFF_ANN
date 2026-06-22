const fs = require("fs");
const path = require("path");
const os = require("os");
const childProcess = require("child_process");

const repo = process.cwd();
const template = path.join(repo, "examples", "classification", "Iris_classification.mlx");

function esc(text) {
  return text.replaceAll("]]>", "]]]]><![CDATA[>");
}

function makeMlx(outPath, code) {
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-safe-"));
  childProcess.execFileSync("unzip", ["-q", template, "-d", work]);
  fs.writeFileSync(
    path.join(work, "matlab", "document.xml"),
    `<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:p><w:pPr><w:pStyle w:val="code"/></w:pPr><w:r><w:t><![CDATA[${esc(code)}]]></w:t></w:r></w:p></w:body></w:document>`
  );
  fs.writeFileSync(
    path.join(work, "matlab", "output.xml"),
    `<?xml version="1.0" encoding="UTF-8"?><embeddedOutputs><metaData><evaluationState>manual</evaluationState><layoutState>code</layoutState><outputStatus>ready</outputStatus></metaData><outputArray type="array"/><regionArray type="array"/></embeddedOutputs>`
  );
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.rmSync(outPath, { force: true });
  childProcess.execFileSync("zip", ["-qr", outPath, "."], { cwd: work });
  fs.rmSync(work, { recursive: true, force: true });
}

const common = `
function repoRoot = locateRepoRoot()
    starts = string.empty;
    scriptPath = mfilename('fullpath');
    if ~isempty(scriptPath)
        starts(end+1) = string(fileparts(scriptPath)); %#ok<AGROW>
    end
    stack = dbstack('-completenames');
    for k = 1:numel(stack)
        if isfield(stack(k), 'file') && ~isempty(stack(k).file)
            starts(end+1) = string(fileparts(stack(k).file)); %#ok<AGROW>
        end
    end
    starts(end+1) = string(pwd);
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
    error('Could not locate release root. Run from inside the release folder.');
end

function addReleasePaths(repoRoot)
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(genpath(fullfile(repoRoot, 'examples')));
    addpath(genpath(fullfile(repoRoot, 'src', 'dynamical_systems')));
end

function assertDisjointSplits(trainIdx, valIdx, testIdx, n, label)
    trainIdx = trainIdx(:); valIdx = valIdx(:); testIdx = testIdx(:);
    assert(~isempty(trainIdx) && ~isempty(valIdx) && ~isempty(testIdx), '%s split is empty.', label);
    allIdx = [trainIdx; valIdx; testIdx];
    assert(all(allIdx >= 1 & allIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(isempty(intersect(trainIdx, testIdx)), '%s train/test leakage.', label);
    assert(isempty(intersect(valIdx, testIdx)), '%s validation/test leakage.', label);
    assert(numel(unique(allIdx)) == numel(allIdx), '%s duplicate split indices.', label);
end

function stats = fitFeatureStandardization(XTrain)
    stats.mu = mean(XTrain, 2);
    stats.sigma = std(XTrain, 0, 2);
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1 / sqrt(size(XTrain, 1));
end

function X = applyFeatureStandardization(X, stats)
    X = stats.scale * ((X - stats.mu) ./ stats.sigma);
    X(~isfinite(X)) = 0;
end

function assertFeatureStatsFromTraining(stats, X, trainIdx, label)
    expected = fitFeatureStandardization(X(:, trainIdx));
    assertStatsClose(stats.mu, expected.mu, label, 'feature mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'feature standard deviation');
    assertStatsClose(stats.scale, expected.scale, label, 'feature scale');
end

function stats = fitImageStandardization(XTrain)
    XTrain = im2single(ensure4d(XTrain));
    stats.mu = mean(XTrain, 4);
    stats.sigma = std(XTrain, 0, 4);
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1 / sqrt(28*28);
end

function X = applyImageStandardization(X, stats)
    X = im2single(ensure4d(X));
    X = stats.scale * ((X - stats.mu) ./ stats.sigma);
    X(~isfinite(X)) = 0;
end

function assertImageStatsFromTraining(stats, X, trainIdx, label)
    expected = fitImageStandardization(X(:, :, :, trainIdx));
    assertStatsClose(stats.mu, expected.mu, label, 'image mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'image standard deviation');
    assertStatsClose(stats.scale, expected.scale, label, 'image scale');
end

function assertStatsClose(actual, expected, label, statName)
    actual = double(actual); expected = double(expected);
    assert(isequal(size(actual), size(expected)), '%s %s size mismatch.', label, statName);
    tol = 1e-10 * max(1, max(abs(expected(:))));
    assert(all(abs(actual(:) - expected(:)) <= tol), '%s %s was not fitted from the training split.', label, statName);
end

function net = loadSavedNetwork(modelPath)
    assert(isfile(modelPath), 'Missing saved network: %s', modelPath);
    S = load(modelPath);
    names = fieldnames(S);
    for candidate = ["net", "dlnet", "trainedNet", "netObj"]
        name = char(candidate);
        if isfield(S, name) && isNetworkLike(S.(name))
            net = S.(name);
            return
        end
    end
    for k = 1:numel(names)
        if isNetworkLike(S.(names{k}))
            net = S.(names{k});
            return
        end
    end
    error('No network-like variable found in %s.', modelPath);
end

function tf = isNetworkLike(v)
    c = lower(class(v));
    tf = contains(c, 'network') || strcmp(c, 'dlnetwork') || strcmp(c, 'seriesnetwork') || strcmp(c, 'dagnetwork');
end

function A = numericOutput(Y)
    if isa(Y, 'dlarray'), A = extractdata(Y); else, A = Y; end
    try, A = gather(A); catch, end
    A = double(A);
end

function Y = predictCB(net, X)
    try
        Y = predict(net, dlarray(single(X), 'CB'));
    catch
        Y = predict(net, single(X.'));
    end
end

function X = ensure4d(X)
    sz = size(X);
    if numel(sz) == 3 && sz(3) ~= 1
        X = reshape(X, 28, 28, 1, sz(3));
    elseif numel(sz) == 2 && sz(1) == 28*28
        X = reshape(X, 28, 28, 1, sz(2));
    elseif numel(sz) == 4
    else
        error('Unsupported image array shape.');
    end
end

function assertNoSharedImages(trainImages, testImages, label)
    trainImages = ensure4d(trainImages);
    testImages = ensure4d(testImages);
    trainRows = reshape(trainImages, [], size(trainImages, 4)).';
    testRows = reshape(testImages, [], size(testImages, 4)).';
    assert(isempty(intersect(trainRows, testRows, 'rows')), '%s exact duplicate images across train/test.', label);
end
`;

const classificationDataHelpers = `
function [X, labels, batchSize] = classificationDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'iris'
            load fisheriris
            labels = categorical(species);
            X = meas.';
            batchSize = 150;
        case 'breast_cancer'
            S = load(fullfile(dataDir, 'breast_cancer_dataset.mat'));
            [labels, X] = breastCancerData(S.data);
            batchSize = 256;
        case 'car_quality'
            S = load(fullfile(dataDir, 'car_dataset.mat'));
            D = S.data; D(~isfinite(D)) = 0;
            labels = categorical(D(:, 7));
            X = D(:, 1:6).';
            batchSize = 256;
        case 'mushroom'
            S = load(fullfile(dataDir, 'mushroom_dataset.mat'));
            D = S.data; D(~isfinite(D)) = 0;
            labels = categorical(D(:, 1));
            X = D(:, 2:end).';
            batchSize = 4096;
        otherwise
            error('Unknown classification task: %s', task);
    end
    X(~isfinite(X)) = 0;
end

function [labels, X] = breastCancerData(data)
    if istable(data)
        labels = categorical(string(data{:, 2})); X = data{:, 3:end}.';
    elseif iscell(data)
        labels = categorical(string(data(:, 2))); X = cellfun(@double, data(:, 3:end)).';
    elseif isnumeric(data)
        [~, ~, g] = unique(data(:, 2)); labels = categorical(g); X = data(:, 3:end).';
    elseif isstruct(data) && isfield(data, 'features') && isfield(data, 'labels')
        X = data.features.'; labels = categorical(data.labels(:));
    else
        error('Unsupported breast-cancer data format.');
    end
    X(~isfinite(X)) = 0;
end
`;

function classificationTrain(task, modelFile, title) {
  return `%% ${title}
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
TASK = '${task}';
MODEL_FILE = '${modelFile}';

clearvars -except repoRoot WIDTH MAX_EPOCHS LR TASK MODEL_FILE
close all; clc; rng(42, 'twister');

[X, labels, batchSize] = classificationDataset(repoRoot, TASK);
[trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), 42, TASK);
featureStats = fitFeatureStandardization(X(:, trainIdx));
assertFeatureStatsFromTraining(featureStats, X, trainIdx, TASK);

XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = labels(trainIdx);
XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = labels(valIdx);
XTest = applyFeatureStandardization(X(:, testIdx), featureStats); YTest = labels(testIdx);
K = numel(categories(YTrain));

layers = [
    featureInputLayer(size(XTrain, 1), Normalization='none')
    fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(K, 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    softmaxLayer];
net = dlnetwork(layers);

options = trainingOptions('adam', ...
    'MaxEpochs', MAX_EPOCHS, 'Metrics', ["accuracy"], 'MiniBatchSize', batchSize, ...
    'InitialLearnRate', LR, 'Shuffle', 'every-epoch', ...
    'ValidationData', {dlarray(single(XVal), 'CB'), YVal}, ...
    'ValidationFrequency', 10, 'Plots', 'training-progress', 'ExecutionEnvironment', 'auto', 'Verbose', 0);

[net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrain, net, 'crossentropy', options);

scores = numericOutput(predictCB(net, XTest));
[~, idx] = max(scores, [], 1);
predictedLabels = categorical(idx(:), 1:K, categories(YTrain));
acc = mean(predictedLabels(:) == YTest(:));
fprintf('%s test accuracy: %.2f%%\\n', TASK, 100*acc);

figure('Color', 'w');
confusionchart(YTest, predictedLabels);
title(sprintf('%s confusion matrix', TASK));

split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
metadata.task = TASK; metadata.width = WIDTH; metadata.normalization = 'train-fitted feature z-score';
outDir = fullfile(repoRoot, 'trained_networks'); if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, MODEL_FILE), 'net', 'info', 'split', 'featureStats', 'metadata', '-v7.3');

${common}

function [trainIdx, valIdx, testIdx] = fixedSplit(n, seed, label)
    rng(seed, 'twister');
    cv = cvpartition(n, 'HoldOut', 0.2);
    trainFullIdx = find(training(cv));
    testIdx = find(test(cv));
    cvVal = cvpartition(numel(trainFullIdx), 'HoldOut', 0.2);
    trainIdx = trainFullIdx(training(cvVal));
    valIdx = trainFullIdx(test(cvVal));
    assertDisjointSplits(trainIdx, valIdx, testIdx, n, label);
end

${classificationDataHelpers}
`;
}

const regressionDataHelpers = `
function [X, y, batchSize] = regressionDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'abalone'
            S = load(fullfile(dataDir, 'abalone_dataset.mat'));
            D = S.data; D(~isfinite(D)) = 0;
            y = D(:, end);
            X = D(:, 1:end-1).';
            batchSize = size(X, 2);
        case 'toyota'
            S = load(fullfile(dataDir, 'toyota_dataset.mat'));
            D = S.data; D(~isfinite(D)) = 0;
            y = D(:, 3);
            X = D(:, [1, 2, 4:size(D, 2)]).';
            batchSize = min(4096, size(X, 2));
        otherwise
            error('Unknown regression task: %s', task);
    end
    X(~isfinite(X)) = 0;
end
`;

function regressionTrain(task, modelFile, splitSeed, title) {
  return `%% ${title}
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
TASK = '${task}';
MODEL_FILE = '${modelFile}';
SPLIT_SEED = ${splitSeed};

clearvars -except repoRoot WIDTH MAX_EPOCHS LR TASK MODEL_FILE SPLIT_SEED
close all; clc; rng(SPLIT_SEED, 'twister');

[X, y, batchSize] = regressionDataset(repoRoot, TASK);
[trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), SPLIT_SEED, TASK);
featureStats = fitFeatureStandardization(X(:, trainIdx));
assertFeatureStatsFromTraining(featureStats, X, trainIdx, TASK);

XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = y(trainIdx);
XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = y(valIdx);
XTest = applyFeatureStandardization(X(:, testIdx), featureStats); YTest = y(testIdx);
targetStats.yMu = mean(YTrain);
targetStats.yStd = std(YTrain); if targetStats.yStd == 0, targetStats.yStd = 1; end
YTrainNorm = (YTrain - targetStats.yMu) / targetStats.yStd;
YValNorm = (YVal - targetStats.yMu) / targetStats.yStd;

layers = [
    featureInputLayer(size(XTrain, 1), Normalization='none')
    fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(1, 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')];
net = dlnetwork(layers);

options = trainingOptions('adam', ...
    'MaxEpochs', MAX_EPOCHS, 'Metrics', "rmse", 'MiniBatchSize', batchSize, ...
    'InitialLearnRate', LR, 'Shuffle', 'every-epoch', ...
    'ValidationData', {dlarray(single(XVal), 'CB'), YValNorm}, ...
    'ValidationFrequency', 10, 'Plots', 'training-progress', 'ExecutionEnvironment', 'auto', 'Verbose', 0);

[net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrainNorm, net, 'mse', options);

YpredNorm = numericOutput(predictCB(net, XTest));
Ypred = YpredNorm(:) * targetStats.yStd + targetStats.yMu;
err = Ypred(:) - double(YTest(:));
fprintf('%s test RMSE: %.4f\\n', TASK, sqrt(mean(err.^2, 'omitnan')));
fprintf('%s signed error mean: %.4f, SD: %.4f\\n', TASK, mean(err, 'omitnan'), std(err, 0, 'omitnan'));

figure('Color', 'w');
histogram(err, 40); grid on; xlabel('Signed error'); ylabel('Count');
title(sprintf('%s signed error distribution', TASK));

split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
metadata.task = TASK; metadata.width = WIDTH; metadata.normalization = 'train-fitted feature z-score';
outDir = fullfile(repoRoot, 'trained_networks'); if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, MODEL_FILE), 'net', 'info', 'split', 'featureStats', 'targetStats', 'metadata', '-v7.3');

${common}

function [trainIdx, valIdx, testIdx] = fixedSplit(n, seed, label)
    rng(seed, 'twister');
    cv = cvpartition(n, 'HoldOut', 0.2);
    trainFullIdx = find(training(cv));
    testIdx = find(test(cv));
    cvVal = cvpartition(numel(trainFullIdx), 'HoldOut', 0.2);
    trainIdx = trainFullIdx(training(cvVal));
    valIdx = trainFullIdx(test(cvVal));
    assertDisjointSplits(trainIdx, valIdx, testIdx, n, label);
end

${regressionDataHelpers}
`;
}

const mnistVariantsTrain = `%% Train MNIST-family classifiers without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

WIDTH = 16000;
NUM_EPOCHS = 1000;
BATCH_SIZE = 500;
LEARN_RATE = 0.01;
NETWORK_DIR = fullfile(repoRoot, 'trained_networks');
DATA_DIR = fullfile(repoRoot, 'data');

if ~exist('target_list', 'var')
    target_list = {'afro_mnist_ethiopic.mat','afro_mnist_nko.mat', ...
                   'afro_mnist_osmanya.mat','afro_mnist_vai.mat', ...
                   'kmnist.mat','mnist_fashion.mat','mnist.mat'};
end

clearvars -except WIDTH NUM_EPOCHS BATCH_SIZE LEARN_RATE NETWORK_DIR DATA_DIR target_list repoRoot
close all; clc; rng(42, 'twister');
if ~exist(NETWORK_DIR, 'dir'), mkdir(NETWORK_DIR); end

for ii = 1:numel(target_list)
    ds = target_list{ii};
    dataPath = fullfile(DATA_DIR, ds);
    assert(isfile(dataPath), 'Missing MNIST-family dataset: %s', dataPath);
    D = load(dataPath);
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, ds);
    [trainIdx, valIdx] = mnistTrainValSplit(numel(YTrainFull), 42, ds);
    imageStats = fitImageStandardization(XTrainFull(:, :, :, trainIdx));
    assertImageStatsFromTraining(imageStats, XTrainFull, trainIdx, ds);

    XTrain = applyImageStandardization(XTrainFull(:, :, :, trainIdx), imageStats); YTrain = YTrainFull(trainIdx);
    XVal = applyImageStandardization(XTrainFull(:, :, :, valIdx), imageStats); YVal = YTrainFull(valIdx);
    XTestEval = applyImageStandardization(XTest, imageStats);

    layers = [
        imageInputLayer([28 28 1], Normalization='none')
        fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
        tanhLayer
        fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
        tanhLayer
        fullyConnectedLayer(numOut, 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
        softmaxLayer];
    net = dlnetwork(layers);
    options = trainingOptions('adam', ...
        'MaxEpochs', NUM_EPOCHS, 'Metrics', ["accuracy"], 'MiniBatchSize', BATCH_SIZE, ...
        'InitialLearnRate', LEARN_RATE, 'Shuffle', 'every-epoch', ...
        'ValidationData', {XVal, YVal}, 'ValidationFrequency', 10, ...
        'Plots', 'training-progress', 'ExecutionEnvironment', 'auto', 'Verbose', false);
    [net, info] = trainnet(XTrain, YTrain, net, 'crossentropy', options);
    scores = numericOutput(predictImageBatches(net, XTestEval, 256));
    predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
    yHat = categorical(predIdx(:), 1:numOut);
    acc = 100 * mean(yHat(:) == YTest(:));
    fprintf('%s test accuracy: %.2f%%\\n', ds, acc);

    split.trainIdx = trainIdx; split.valIdx = valIdx; split.testSplit = 'official_test_struct';
    metadata.task = erase(ds, '.mat'); metadata.dataset = ds; metadata.width = WIDTH; metadata.normalization = 'train-fitted image z-score';
    save(fullfile(NETWORK_DIR, sprintf('%s_network.mat', metadata.task)), 'net', 'info', 'acc', 'split', 'imageStats', 'metadata', '-v7.3');
end

${common}

function [XTrain, YTrain, XTest, YTest, numOut] = mnistData(training, test, label)
    assert(isfield(training, 'images') && isfield(training, 'labels'), '%s missing training fields.', label);
    assert(isfield(test, 'images') && isfield(test, 'labels'), '%s missing test fields.', label);
    assertNoSharedImages(training.images, test.images, label);
    XTrain = ensure4d(training.images);
    XTest = ensure4d(test.images);
    YTrain = categorical(double(training.labels(:)) + 1);
    YTest = categorical(double(test.labels(:)) + 1);
    numOut = max(double(training.labels)) + 1;
end

function [trainIdx, valIdx] = mnistTrainValSplit(n, seed, label)
    rng(seed, 'twister');
    idx = randperm(n);
    nTrain = round(0.8 * n);
    trainIdx = idx(1:nTrain);
    valIdx = idx(nTrain+1:end);
    assert(~isempty(trainIdx) && ~isempty(valIdx), '%s train/validation split is empty.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
end

function Y = predictImageBatches(net, X, batchSize)
    nObs = size(X, 4);
    Y = [];
    for first = 1:batchSize:nObs
        last = min(first + batchSize - 1, nObs);
        YBatch = numericOutput(predict(net, X(:, :, :, first:last)));
        if isempty(Y)
            if size(YBatch, 1) == numel(first:last)
                Y = zeros(nObs, size(YBatch, 2));
            else
                Y = zeros(size(YBatch, 1), nObs);
            end
        end
        if size(Y, 1) == nObs
            Y(first:last, :) = YBatch;
        else
            Y(:, first:last) = YBatch;
        end
    end
end

function idx = scoresToClassIndex(scores, nClasses, nObs)
    if size(scores, 1) == nObs && size(scores, 2) == nClasses
        [~, idx] = max(scores, [], 2);
    elseif size(scores, 1) == nClasses && size(scores, 2) == nObs
        [~, idx] = max(scores, [], 1); idx = idx(:);
    else
        error('Unexpected score size %s.', mat2str(size(scores)));
    end
end
`;

function mnistSingleTrain(dataset) {
  return mnistVariantsTrain
    .replace("%% Train MNIST-family classifiers without normalization leakage", `%% Train ${dataset.replace(".mat", "")} classifier without normalization leakage`)
    .replace(`if ~exist('target_list', 'var')
    target_list = {'afro_mnist_ethiopic.mat','afro_mnist_nko.mat', ...
                   'afro_mnist_osmanya.mat','afro_mnist_vai.mat', ...
                   'kmnist.mat','mnist_fashion.mat','mnist.mat'};
end`, `target_list = {'${dataset}'};`);
}

function classificationTest() {
  return `%% Test saved tabular classification networks without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

tasks = struct( ...
    'task', {'iris', 'breast_cancer', 'car_quality', 'mushroom'}, ...
    'label', {'Iris', 'Breast cancer', 'Car quality', 'Mushroom'}, ...
    'model', {'Iris_network.mat', 'BC_network.mat', 'Car_quality_network.mat', 'Mushroom_network.mat'});

for t = 1:numel(tasks)
    task = tasks(t);
    modelPath = fullfile(repoRoot, 'trained_networks', task.model);
    Snet = load(modelPath);
    assert(isfield(Snet, 'split'), '%s lacks saved split metadata; rerun the training script.', task.model);
    assert(isfield(Snet, 'featureStats'), '%s lacks train-fitted featureStats; rerun the training script.', task.model);
    net = loadSavedNetwork(modelPath);
    [X, labels, ~] = classificationDataset(repoRoot, task.task);
    assertDisjointSplits(Snet.split.trainIdx, Snet.split.valIdx, Snet.split.testIdx, size(X, 2), task.label);
    assertFeatureStatsFromTraining(Snet.featureStats, X, Snet.split.trainIdx, task.label);
    XTest = applyFeatureStandardization(X(:, Snet.split.testIdx), Snet.featureStats);
    YTest = labels(Snet.split.testIdx);
    YTrain = labels(Snet.split.trainIdx);
    scores = numericOutput(predictCB(net, XTest));
    [~, idx] = max(scores, [], 1);
    yHat = categorical(idx(:), 1:numel(categories(YTrain)), categories(YTrain));
    fprintf('%s saved-network accuracy: %.2f%%\\n', task.label, 100*mean(yHat(:) == YTest(:)));
end

${common}
${classificationDataHelpers}
`;
}

function regressionTest() {
  return `%% Test saved regression networks without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

tasks = struct( ...
    'task', {'abalone', 'toyota'}, ...
    'label', {'Abalone age', 'Toyota price'}, ...
    'model', {'Abalone_network.mat', 'Toyota_network.mat'});

for t = 1:numel(tasks)
    task = tasks(t);
    modelPath = fullfile(repoRoot, 'trained_networks', task.model);
    Snet = load(modelPath);
    assert(isfield(Snet, 'split'), '%s lacks saved split metadata; rerun the training script.', task.model);
    assert(isfield(Snet, 'featureStats'), '%s lacks train-fitted featureStats; rerun the training script.', task.model);
    assert(isfield(Snet, 'targetStats'), '%s lacks targetStats; rerun the training script.', task.model);
    net = loadSavedNetwork(modelPath);
    [X, y, ~] = regressionDataset(repoRoot, task.task);
    assertDisjointSplits(Snet.split.trainIdx, Snet.split.valIdx, Snet.split.testIdx, size(X, 2), task.label);
    assertFeatureStatsFromTraining(Snet.featureStats, X, Snet.split.trainIdx, task.label);
    XTest = applyFeatureStandardization(X(:, Snet.split.testIdx), Snet.featureStats);
    Ytrue = double(y(Snet.split.testIdx));
    YpredNorm = numericOutput(predictCB(net, XTest));
    Ypred = YpredNorm(:) * Snet.targetStats.yStd + Snet.targetStats.yMu;
    err = Ypred(:) - Ytrue(:);
    [r, p] = corrStats(Ypred(:), Ytrue(:));
    fprintf('%s signed error mean: %.4f\\n', task.label, mean(err, 'omitnan'));
    fprintf('%s signed error SD: %.4f\\n', task.label, std(err, 0, 'omitnan'));
    fprintf('%s RMSE: %.4f\\n', task.label, sqrt(mean(err.^2, 'omitnan')));
    fprintf('%s Pearson r: %.4f, p: %.4g\\n', task.label, r, p);
    figure('Color', 'w'); histogram(err, 40); grid on; xlabel('Signed error'); ylabel('Count'); title(sprintf('%s signed error distribution', task.label));
end

${common}
${regressionDataHelpers}

function [r, p] = corrStats(yPred, yTrue)
    valid = isfinite(yPred) & isfinite(yTrue);
    yPred = yPred(valid); yTrue = yTrue(valid);
    if numel(yPred) < 3 || std(yPred) == 0 || std(yTrue) == 0
        r = NaN; p = NaN; return
    end
    [R, P] = corrcoef(yPred, yTrue); r = R(1, 2); p = P(1, 2);
end
`;
}

function mnistSavedTest() {
  return `%% Test saved MNIST-family networks without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

networkFiles = [
    "mnist_network.mat"
    "kmnist_network.mat"
    "mnist_fashion_network.mat"
    "afro_mnist_ethiopic_network.mat"
    "afro_mnist_nko_network.mat"
    "afro_mnist_osmanya_network.mat"
    "afro_mnist_vai_network.mat"
];

for i = 1:numel(networkFiles)
    networkPath = fullfile(repoRoot, 'trained_networks', networkFiles(i));
    Snet = load(networkPath);
    assert(isfield(Snet, 'split') && strcmp(Snet.split.testSplit, 'official_test_struct'), '%s lacks official test-split metadata; rerun MNIST_variants.', networkFiles(i));
    assert(isfield(Snet, 'imageStats'), '%s lacks train-fitted imageStats; rerun MNIST_variants.', networkFiles(i));
    net = loadSavedNetwork(networkPath);
    dataFile = erase(networkFiles(i), '_network.mat') + ".mat";
    D = load(fullfile(repoRoot, 'data', dataFile));
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, dataFile);
    assert(isfield(Snet.split, 'trainIdx') && isfield(Snet.split, 'valIdx'), '%s lacks train/validation split metadata.', networkFiles(i));
    assert(isempty(intersect(Snet.split.trainIdx(:), Snet.split.valIdx(:))), '%s train/validation leakage.', networkFiles(i));
    assert(all(Snet.split.trainIdx(:) >= 1 & Snet.split.trainIdx(:) <= numel(YTrainFull)), '%s train indices out of range.', networkFiles(i));
    assert(all(Snet.split.valIdx(:) >= 1 & Snet.split.valIdx(:) <= numel(YTrainFull)), '%s validation indices out of range.', networkFiles(i));
    assertImageStatsFromTraining(Snet.imageStats, XTrainFull, Snet.split.trainIdx, dataFile);
    XTestEval = applyImageStandardization(XTest, Snet.imageStats);
    scores = numericOutput(predictImageBatches(net, XTestEval, 256));
    predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
    yHat = categorical(predIdx(:), 1:numOut);
    fprintf('%s test-split accuracy: %.2f%%\\n', networkFiles(i), 100*mean(yHat(:) == YTest(:)));
end

${common}

function [XTrain, YTrain, XTest, YTest, numOut] = mnistData(training, test, label)
    assertNoSharedImages(training.images, test.images, label);
    XTrain = ensure4d(training.images);
    XTest = ensure4d(test.images);
    YTrain = categorical(double(training.labels(:)) + 1);
    YTest = categorical(double(test.labels(:)) + 1);
    numOut = max(double(training.labels)) + 1;
end

function Y = predictImageBatches(net, X, batchSize)
    nObs = size(X, 4);
    Y = [];
    for first = 1:batchSize:nObs
        last = min(first + batchSize - 1, nObs);
        YBatch = numericOutput(predict(net, X(:, :, :, first:last)));
        if isempty(Y)
            if size(YBatch, 1) == numel(first:last), Y = zeros(nObs, size(YBatch, 2)); else, Y = zeros(size(YBatch, 1), nObs); end
        end
        if size(Y, 1) == nObs, Y(first:last, :) = YBatch; else, Y(:, first:last) = YBatch; end
    end
end

function idx = scoresToClassIndex(scores, nClasses, nObs)
    if size(scores, 1) == nObs && size(scores, 2) == nClasses
        [~, idx] = max(scores, [], 2);
    elseif size(scores, 1) == nClasses && size(scores, 2) == nObs
        [~, idx] = max(scores, [], 1); idx = idx(:);
    else
        error('Unexpected score size %s.', mat2str(size(scores)));
    end
end
`;
}

makeMlx(path.join(repo, "examples", "classification", "Iris_classification.mlx"), classificationTrain("iris", "Iris_network.mat", "Train Iris classifier without normalization leakage"));
makeMlx(path.join(repo, "examples", "classification", "Breast_cancer_classification.mlx"), classificationTrain("breast_cancer", "BC_network.mat", "Train breast-cancer classifier without normalization leakage"));
makeMlx(path.join(repo, "examples", "classification", "Car_quality_classification.mlx"), classificationTrain("car_quality", "Car_quality_network.mat", "Train car-quality classifier without normalization leakage"));
makeMlx(path.join(repo, "examples", "classification", "Mushroom_classification.mlx"), classificationTrain("mushroom", "Mushroom_network.mat", "Train mushroom classifier without normalization leakage"));
makeMlx(path.join(repo, "examples", "regression", "Abalone_age_prediction.mlx"), regressionTrain("abalone", "Abalone_network.mat", 1, "Train abalone age regressor without normalization leakage"));
makeMlx(path.join(repo, "examples", "regression", "Car_price_regression.mlx"), regressionTrain("toyota", "Toyota_network.mat", 42, "Train Toyota price regressor without normalization leakage"));
makeMlx(path.join(repo, "examples", "classification", "MNIST_variants.mlx"), mnistVariantsTrain);
for (const dataset of [
  "afro_mnist_ethiopic.mat",
  "afro_mnist_nko.mat",
  "afro_mnist_osmanya.mat",
  "afro_mnist_vai.mat",
  "kmnist.mat",
  "mnist_fashion.mat",
  "mnist.mat",
]) {
  const task = dataset.replace(".mat", "");
  makeMlx(path.join(repo, "examples", "classification", `MNIST_${task}.mlx`), mnistSingleTrain(dataset));
}

makeMlx(path.join(repo, "tests", "test_iris_saved_network.mlx"), classificationTest());
makeMlx(path.join(repo, "tests", "test_tabular_saved_networks.mlx"), classificationTest());
makeMlx(path.join(repo, "tests", "test_regression_saved_networks.mlx"), regressionTest());
makeMlx(path.join(repo, "tests", "test_mnist_saved_networks.mlx"), mnistSavedTest());

const dynamicalHelpers = `
function stats = fitStateStandardization(xTrain)
    stats.mu = mean(xTrain, 1);
    stats.sigma = std(xTrain, 0, 1);
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1 / sqrt(size(xTrain, 2));
end

function x = applyStateStandardization(x, stats)
    x = stats.scale * ((x - stats.mu) ./ stats.sigma);
    x(~isfinite(x)) = 0;
end

function assertStateStatsFromTraining(stats, xTrain, label)
    expected = fitStateStandardization(xTrain);
    assertStatsClose(stats.mu, expected.mu, label, 'state mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'state standard deviation');
    assertStatsClose(stats.scale, expected.scale, label, 'state scale');
end

function yhat = predictDynamicalStep(net, inCol)
    try
        yhat = predict(net, inCol(:).');
    catch
        yhat = predictCB(net, inCol(:));
    end
    yhat = numericOutput(yhat);
end
`;

const dynamicalTrain = `%% Train dynamical-system NAR network without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

ii = 5;
perturbationScale = 0.01;

clearvars -except repoRoot ii perturbationScale
close all; clc;
dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
dynamicsName = string(dyn{ii, 1});
ic = dyn{ii, 4:6}.';
d = dyn{ii, 2};
tSpan = [0, 3000];
inputDelay = 1;

[tTrain, xRaw] = ode45(@(t,y) int_dyn(y, dynamicsName, 0, 0, 0, 'Simulate'), tSpan, ic(1:d));
[inpRaw, tgtRaw] = extractTimeDelayedInputs(xRaw, inputDelay);
N = size(inpRaw, 3);
idx = floor(0.8 * N);
assert(idx > 0 && idx < N, 'Dynamical-system train/validation split is empty.');
trainStateRows = 1:(idx + inputDelay);
valStateRows = (idx + 1 + inputDelay):size(xRaw, 1);
assert(isempty(intersect(trainStateRows(:), valStateRows(:))), 'Dynamical-system train/validation leakage.');
stateStats = fitStateStandardization(xRaw(trainStateRows, :));
assertStateStatsFromTraining(stateStats, xRaw(trainStateRows, :), dynamicsName);
xTrue = applyStateStandardization(xRaw, stateStats);

figure('Color', 'w');
subplot(2,1,1); plot(tTrain, xTrue, 'LineWidth', 1); title('True Time Series'); xlabel('Time'); ylabel('State');
if size(xTrue, 2) > 1
    subplot(2,1,2); plot(xTrue(:,1), xTrue(:,2), 'LineWidth', 1); title('Phase Plot'); xlabel('x_1'); ylabel('x_2');
end

[inp, tgt] = extractTimeDelayedInputs(xTrue, inputDelay);
X = dlarray(inp, 'CTB'); Y = dlarray(tgt, 'CTB');
XTrain = X(:, :, 1:idx); YTrain = Y(:, :, 1:idx);
XVal = X(:, :, idx+1:end); YVal = Y(:, :, idx+1:end);

layers = [
    sequenceInputLayer(size(XTrain, 1))
    fullyConnectedLayer(16000, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(16000, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
    tanhLayer
    fullyConnectedLayer(size(xTrue, 2), 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')];
net = dlnetwork(layerGraph(layers));

options = trainingOptions('adam', ...
    MaxEpochs=20000, InitialLearnRate=0.01, MiniBatchSize=size(XTrain, 3), ...
    Shuffle='never', Plots='training-progress', ValidationData={XVal, YVal}, ...
    ValidationFrequency=10, Verbose=0, L2Regularization=0, ...
    GradientDecayFactor=0.99, SquaredGradientDecayFactor=0.999, Epsilon=5e-6);
[net, trainingInfo] = trainnet(XTrain, YTrain, net, 'mse', options);

icTest = ic + perturbationScale * randn(size(ic));
assert(norm(icTest - ic) > 0, 'Dynamical-system test initial condition was reused.');
[~, xTestRaw] = ode45(@(t,y) int_dyn(y, dynamicsName, 0, 0, 0, 'Simulate'), 5*tSpan, icTest(1:d));
xTestTrue = applyStateStandardization(xTestRaw, stateStats);
xPred = xTestTrue(1:inputDelay, :);
for k = 1:(size(xTestTrue, 1) - inputDelay)
    Xf = extractTimeDelayedInputsFeedback(xPred, inputDelay);
    yP = predictDynamicalStep(net, Xf(:, end));
    xPred = [xPred; yP(:).']; %#ok<AGROW>
end

figure('Color', 'w');
if size(xTrue, 2) > 1
    plot(xTestTrue(:,1), xTestTrue(:,2), 'k', 'LineWidth', 1); hold on;
    plot(xPred(:,1), xPred(:,2), 'r--', 'LineWidth', 1); hold off;
    legend('True', 'Predicted'); title('Leakage-safe perturbed rollout'); xlabel('x_1'); ylabel('x_2');
end

split.trainTimeRows = trainStateRows; split.valTimeRows = valStateRows; split.inputDelay = inputDelay;
metadata.task = char(dynamicsName); metadata.normalization = 'train-segment-fitted state z-score';
outDir = fullfile(repoRoot, 'trained_networks'); if ~exist(outDir, 'dir'), mkdir(outDir); end
save(fullfile(outDir, sprintf('%s_network.mat', dynamicsName)), 'net', 'trainingInfo', 'inputDelay', 'split', 'stateStats', 'metadata', '-v7.3');

${common}
${dynamicalHelpers}
`;

function dynamicalSingleTrain(row, name) {
  return dynamicalTrain
    .replace("%% Train dynamical-system NAR network without normalization leakage", `%% Train ${name} dynamical-system NAR network without normalization leakage`)
    .replace("ii = 5;", `ii = ${row};`);
}

const dynamicalTest = `%% Test saved dynamical-system networks without normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

modelFiles = [
    "Lorenz_network.mat"
    "MO0_network.mat"
    "MO5_network.mat"
    "MO7_network.mat"
    "MO13_network.mat"
    "Rikitake_network.mat"
    "SprottB_network.mat"
    "SprottC_network.mat"
    "SprottS_network.mat"
];

for i = 1:numel(modelFiles)
    modelPath = fullfile(repoRoot, 'trained_networks', modelFiles(i));
    Snet = load(modelPath);
    assert(isfield(Snet, 'stateStats'), '%s lacks train-fitted stateStats; rerun the dynamical-system training script.', modelFiles(i));
    assert(isfield(Snet, 'split'), '%s lacks split metadata; rerun the dynamical-system training script.', modelFiles(i));
    net = loadSavedNetwork(modelPath);
    dynamicsName = erase(modelFiles(i), '_network.mat');
    rmse = dynamicalSystemTestRmse(repoRoot, net, Snet, dynamicsName);
    fprintf('%s perturbed test trajectory RMSE: %.4g\\n', dynamicsName, rmse);
end

${common}
${dynamicalHelpers}

function rmse = dynamicalSystemTestRmse(repoRoot, net, Snet, dynamicsName)
    dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
    names = string(dyn{:, 1});
    row = find(names == string(dynamicsName), 1);
    assert(~isempty(row), 'No dynamics-list row found for %s.', dynamicsName);
    d = dyn{row, 2};
    ic = dyn{row, 4:6}.';
    tSpan = [0, 150];
    [~, xTrainRaw] = ode45(@(t,y) int_dyn(y, string(dynamicsName), 0, 0, 0, 'Simulate'), [0, 3000], ic(1:d));
    assertStateStatsFromTraining(Snet.stateStats, xTrainRaw(Snet.split.trainTimeRows, :), dynamicsName);
    rng(9000 + row, 'twister');
    icTest = ic + 0.01 * randn(size(ic));
    assert(norm(icTest - ic) > 0, '%s perturbation failed.', dynamicsName);
    [~, xTestRaw] = ode45(@(t,y) int_dyn(y, string(dynamicsName), 0, 0, 0, 'Simulate'), tSpan, icTest(1:d));
    xTestTrue = applyStateStandardization(xTestRaw, Snet.stateStats);
    inputDelay = Snet.inputDelay;
    xPred = xTestTrue(1:inputDelay, :);
    for k = 1:(size(xTestTrue, 1) - inputDelay)
        Xf = extractTimeDelayedInputsFeedback(xPred, inputDelay);
        yhat = predictDynamicalStep(net, Xf(:, end));
        xPred = [xPred; yhat(:).']; %#ok<AGROW>
    end
    err = xPred - xTestTrue;
    rmse = sqrt(mean(err(:).^2, 'omitnan'));
end
`;

makeMlx(path.join(repo, "examples", "dynamical_systems", "Dynamical_systems.mlx"), dynamicalTrain);
for (const [row, name] of [
  [5, "Lorenz"],
  [34, "MO0"],
  [39, "MO5"],
  [41, "MO7"],
  [47, "MO13"],
  [32, "Rikitake"],
  [8, "SprottB"],
  [9, "SprottC"],
  [25, "SprottS"],
]) {
  makeMlx(path.join(repo, "examples", "dynamical_systems", `train_${name}.mlx`), dynamicalSingleTrain(row, name));
}
makeMlx(path.join(repo, "examples", "dynamical_systems", "test_dynamical_systems.mlx"), dynamicalTest);
makeMlx(path.join(repo, "tests", "test_dynamical_system_saved_networks.mlx"), dynamicalTest);

const motorControlTrain = `%% Train LQR two-link arm network without episode or normalization leakage
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
NETWORK_PATH = fullfile(repoRoot, 'trained_networks', 'LQR_network.mat');
EPISODE_PATH = fullfile(repoRoot, 'data', 'LQR_SUPERVISED_EPISODES_20251002_172744.mat');

clearvars -except repoRoot NETWORK_PATH EPISODE_PATH; clc; close all; rng(1, 'twister');

[A, B, K, dt, T, t, l1, l2] = lqrArmSetup();
numEpisodes = 250;
sampleInterval = 1;
inputsEpisodes = cell(numEpisodes, 1);
outputsEpisodes = cell(numEpisodes, 1);
targetEE = zeros(numEpisodes, 2);

for ep = 1:numEpisodes
    target = randomReachTarget(l1, l2);
    targetEE(ep, :) = target(:).';
    xTargetJoint = inverseKinematicsTarget(target, l1, l2);
    xSim = zeros(4, numel(t));
    epiInputs = [];
    epiOutputs = [];
    consecutiveCounter = 0;
    for k = 1:numel(t)-1
        eJoint = xSim(:, k) - xTargetJoint;
        u = -K * eJoint;
        xSim(:, k+1) = xSim(:, k) + dt*(A*xSim(:, k) + B*u);
        if mod(k, sampleInterval) == 0
            epiInputs = [epiInputs, encodeJointError(eJoint)]; %#ok<AGROW>
            epiOutputs = [epiOutputs, u]; %#ok<AGROW>
        end
        if endEffectorDistance(xSim(:, k+1), target, l1, l2) < 0.01
            consecutiveCounter = consecutiveCounter + 1;
        else
            consecutiveCounter = 0;
        end
        if consecutiveCounter >= 10, break; end
    end
    inputsEpisodes{ep} = epiInputs;
    outputsEpisodes{ep} = epiOutputs;
end

idxEp = randperm(numEpisodes);
numTrainEp = round(0.8 * numEpisodes);
trainEps = idxEp(1:numTrainEp);
valEps = idxEp(numTrainEp+1:end);
assertEpisodeSplit(trainEps, valEps, numEpisodes, 'LQR');

XTrainAll = localConcat(inputsEpisodes, trainEps);
YTrainAll = localConcat(outputsEpisodes, trainEps);
XValAll = localConcat(inputsEpisodes, valEps);
YValAll = localConcat(outputsEpisodes, valEps);
assert(~isempty(XTrainAll) && ~isempty(XValAll), 'LQR train/validation samples are empty.');

mu_X = mean(XTrainAll, 2);
sigma_X = std(XTrainAll, 0, 2); sigma_X(sigma_X == 0) = 1;
mu_Y = mean(YTrainAll, 2);
sigma_Y = std(YTrainAll, 0, 2); sigma_Y(sigma_Y == 0) = 1;
assertStatsClose(mu_X, mean(XTrainAll, 2), 'LQR', 'input mean');
assertStatsClose(sigma_X, replaceZeroStd(std(XTrainAll, 0, 2)), 'LQR', 'input std');
assertStatsClose(mu_Y, mean(YTrainAll, 2), 'LQR', 'output mean');
assertStatsClose(sigma_Y, replaceZeroStd(std(YTrainAll, 0, 2)), 'LQR', 'output std');

numInputs = size(XTrainAll, 1);
XTrain = ((XTrainAll - mu_X) ./ (sigma_X * sqrt(numInputs))).';
YTrain = ((YTrainAll - mu_Y) ./ sigma_Y).';
XVal = ((XValAll - mu_X) ./ (sigma_X * sqrt(numInputs))).';
YVal = ((YValAll - mu_Y) ./ sigma_Y).';

layers = [
    featureInputLayer(numInputs, 'Normalization', 'none', 'Name', 'input')
    fullyConnectedLayer(16e3, 'Name', 'fc1', 'WeightsInitializer', @customWeights1_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')
    tanhLayer('Name', 'relu1')
    fullyConnectedLayer(16e3, 'Name', 'fc2', 'WeightsInitializer', @customWeights2_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')
    tanhLayer('Name', 'relu2')
    fullyConnectedLayer(2, 'Name', 'fc3', 'WeightsInitializer', @customWeights3_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')];
net = dlnetwork(layers);
options = trainingOptions('adam', ...
    'Verbose', 0, 'MaxEpochs', 1e4, 'Metrics', 'rmse', 'InitialLearnRate', 0.01, ...
    'L2Regularization', 0, 'Plots', 'training-progress', 'MiniBatchSize', 5000, ...
    'ValidationData', {XVal, YVal}, 'ValidationFrequency', 10);
[net, tr] = trainnet(XTrain, YTrain, net, 'mse', options);

numSims = 100;
testTargets = zeros(numSims, 2);
for sim = 1:numSims
    testTargets(sim, :) = randomReachTarget(l1, l2).';
end
assertNoSharedTargets(targetEE(trainEps, :), testTargets, 'LQR train/test target leakage.');

testSequence.dt = dt; testSequence.T = T; testSequence.t = t; testSequence.l1 = l1; testSequence.l2 = l2;
testSequence.A = A; testSequence.B = B; testSequence.K = K;
testSequence.mu_X = mu_X; testSequence.sigma_X = sigma_X; testSequence.mu_Y = mu_Y; testSequence.sigma_Y = sigma_Y;
testSequence.numInputs = numInputs; testSequence.numSims = numSims; testSequence.sampleInterval = sampleInterval;
testSequence.trainTargets = targetEE(trainEps, :); testSequence.testTargets = testTargets;
split.trainEps = trainEps; split.valEps = valEps;
metadata.task = 'LQR_two_link_arm'; metadata.normalization = 'train-episode-fitted input/output z-score';

save(NETWORK_PATH, 'net', 'tr', 'mu_X', 'sigma_X', 'mu_Y', 'sigma_Y', 'split', 'testSequence', 'metadata', '-v7.3');
save(EPISODE_PATH, 'inputsEpisodes', 'outputsEpisodes', 'trainEps', 'valEps', 'targetEE', 'testSequence', '-v7.3');

${common}
${lqrHelpers()}
`;

function lqrHelpers() {
  return `
function [A, B, K, dt, T, t, l1, l2] = lqrArmSetup()
    A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    B = [0 0; 0 0; 1 0; 0 1];
    Q = diag([100, 100, 1, 1]); R = diag([0.1, 0.1]);
    K = lqr(A, B, Q, R);
    dt = 0.01; T = 5; t = 0:dt:T; l1 = 1.0; l2 = 0.8;
end

function target = randomReachTarget(l1, l2)
    r = 0.5 + (l1 + l2 - 0.5) * rand;
    angle = 2*pi*rand;
    target = [r*cos(angle); r*sin(angle)];
end

function xTargetJoint = inverseKinematicsTarget(target, l1, l2)
    x = target(1); y = target(2);
    c2 = (x^2 + y^2 - l1^2 - l2^2) / (2*l1*l2);
    c2 = max(min(c2, 1), -1);
    theta2 = acos(c2);
    k1 = l1 + l2*cos(theta2); k2 = l2*sin(theta2);
    theta1 = atan2(y, x) - atan2(k2, k1);
    xTargetJoint = [theta1; theta2; 0; 0];
end

function phi = encodeJointError(eJoint)
    phi = [sin(eJoint(1)); cos(eJoint(1)); sin(eJoint(2)); cos(eJoint(2)); eJoint(3); eJoint(4)];
end

function dist = endEffectorDistance(xState, target, l1, l2)
    th1 = xState(1); th2 = xState(2);
    endEffector = [l1*cos(th1) + l2*cos(th1 + th2); l1*sin(th1) + l2*sin(th1 + th2)];
    dist = hypot(endEffector(1) - target(1), endEffector(2) - target(2));
end

function M = localConcat(C, idx)
    M = [];
    for ii = 1:numel(idx)
        ci = C{idx(ii)};
        if ~isempty(ci), M = [M, ci]; end %#ok<AGROW>
    end
end

function assertEpisodeSplit(trainEps, valEps, n, label)
    assert(~isempty(trainEps) && ~isempty(valEps), '%s episode split is empty.', label);
    assert(all(trainEps >= 1 & trainEps <= n) && all(valEps >= 1 & valEps <= n), '%s episode indices out of range.', label);
    assert(isempty(intersect(trainEps, valEps)), '%s train/validation episode leakage.', label);
    assert(numel(unique([trainEps(:); valEps(:)])) == n, '%s duplicate or missing episode indices.', label);
end

function v = replaceZeroStd(v)
    v(v == 0) = 1;
end

function assertNoSharedTargets(trainTargets, testTargets, message)
    sharedTargets = intersect(round(trainTargets, 12), round(testTargets, 12), 'rows');
    assert(isempty(sharedTargets), message);
end
`;
}

const pongTrain = `%% Train Pong network with hard train/validation split checks
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
NETWORK_PATH = fullfile(repoRoot, 'trained_networks', 'Pong_network.mat');

clearvars -except repoRoot NETWORK_PATH; clc; rng(1, 'twister');
ballSpeed = 0.02;
maxBounceAngle = pi/4;
numSamples = 1000;
[features, labels] = generateTrainingData(numSamples, ballSpeed);

numTrain = round(0.8 * numSamples);
trainIdx = 1:numTrain;
valIdx = (numTrain+1):numSamples;
assertNoOverlap(trainIdx, valIdx, numSamples, 'Pong');
trainFeatures = features(trainIdx, :); trainLabels = labels(trainIdx);
valFeatures = features(valIdx, :); valLabels = labels(valIdx);

layers = createPongNetwork();
options = trainingOptions('adam', ...
    'MaxEpochs', 20000, 'MiniBatchSize', 1000, 'ValidationData', {valFeatures, valLabels}, ...
    'Plots', 'training-progress', 'Verbose', 0, 'Metrics', 'accuracy');
net = dlnetwork(layers);
[net, tr] = trainnet(trainFeatures, trainLabels, net, 'crossentropy', options);

scores = numericOutput(predict(net, valFeatures));
predictedIdx = scoresToClassIndex(scores, 3, numel(valLabels));
accuracy = mean(predictedIdx(:) == grp2idx(valLabels(:)));
fprintf('Pong validation accuracy: %.2f%%\\n', 100*accuracy);

split.trainIdx = trainIdx; split.valIdx = valIdx;
metadata.task = 'Pong'; metadata.normalization = 'none'; metadata.ballSpeed = ballSpeed; metadata.maxBounceAngle = maxBounceAngle;
save(NETWORK_PATH, 'net', 'tr', 'split', 'metadata', 'ballSpeed', 'maxBounceAngle', '-v7.3');

${common}

function [features, labels] = generateTrainingData(numSamples, ballSpeed)
    features = zeros(numSamples, 5);
    threshold = 0.1; centerPos = 0.5;
    labels = categorical(zeros(numSamples, 1), [0 1 2], ["Up", "Down", "Stay"]);
    for i = 1:numSamples
        ballX = rand; ballY = rand; paddleY = rand;
        theta = 2*pi*rand;
        ballVelX = ballSpeed * cos(theta); ballVelY = ballSpeed * sin(theta);
        features(i, :) = [ballX, ballY, ballVelX, ballVelY, paddleY];
        if ballX < 0.7
            if paddleY < centerPos - threshold, action = "Down";
            elseif paddleY > centerPos + threshold, action = "Up";
            else, action = "Stay"; end
        else
            if ballY < paddleY - threshold/2, action = "Up";
            elseif ballY > paddleY + threshold/2, action = "Down";
            else, action = "Stay"; end
        end
        labels(i) = categorical(action, ["Up", "Down", "Stay"]);
    end
end

function layers = createPongNetwork()
    layers = [
        featureInputLayer(5, "Name", "input")
        fullyConnectedLayer(16e3, "Name", "fc1", 'WeightsInitializer', @customWeights1_git, 'WeightLearnRateFactor', 0)
        tanhLayer("Name", "relu1")
        fullyConnectedLayer(16e3, "Name", "fc2", 'WeightsInitializer', @customWeights2_git, 'WeightLearnRateFactor', 0)
        tanhLayer("Name", "relu2")
        fullyConnectedLayer(3, "Name", "fc3", 'WeightsInitializer', @customWeights3_git, 'WeightLearnRateFactor', 0)
        softmaxLayer("Name", "softmax")];
end

function assertNoOverlap(trainIdx, valIdx, n, label)
    assert(~isempty(trainIdx) && ~isempty(valIdx), '%s train/validation split is empty.', label);
    assert(all(trainIdx >= 1 & trainIdx <= n) && all(valIdx >= 1 & valIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(numel(unique([trainIdx(:); valIdx(:)])) == n, '%s duplicate or missing split indices.', label);
end

function idx = scoresToClassIndex(scores, nClasses, nObs)
    if size(scores, 1) == nObs && size(scores, 2) == nClasses
        [~, idx] = max(scores, [], 2);
    elseif size(scores, 1) == nClasses && size(scores, 2) == nObs
        [~, idx] = max(scores, [], 1);
        idx = idx(:);
    else
        error('Unexpected score size %s.', mat2str(size(scores)));
    end
end
`;

const pongTest = `%% Test saved Pong network with hard metadata checks
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
modelPath = fullfile(repoRoot, 'trained_networks', 'Pong_network.mat');
S = load(modelPath);
assert(isfield(S, 'split') && isfield(S.split, 'trainIdx') && isfield(S.split, 'valIdx'), 'Pong network lacks split metadata; rerun Pong.mlx.');
assert(isfield(S, 'metadata') && strcmp(S.metadata.normalization, 'none'), 'Pong network metadata must record no normalization.');
nObs = max([S.split.trainIdx(:); S.split.valIdx(:)]);
assertNoOverlap(S.split.trainIdx, S.split.valIdx, nObs, 'Pong');
net = loadSavedNetwork(modelPath);
assertNetworkIO(net, 5, 3, 'Pong');
fprintf('Pong saved network metadata and split checks passed.\\n');

${common}

function assertNoOverlap(trainIdx, valIdx, n, label)
    assert(~isempty(trainIdx) && ~isempty(valIdx), '%s train/validation split is empty.', label);
    assert(all(trainIdx >= 1 & trainIdx <= n) && all(valIdx >= 1 & valIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(numel(unique([trainIdx(:); valIdx(:)])) == n, '%s duplicate or missing split indices.', label);
end

function assertNetworkIO(net, expectedInputs, expectedOutputs, label)
    y = numericOutput(predict(net, zeros(1, expectedInputs, 'single')));
    assert(numel(y) == expectedOutputs, '%s output size mismatch.', label);
end
`;

const pongSimulationTest = `%% Test saved Pong network against a mirroring opponent
% The left paddle mirrors the ball position. The saved network controls the
% right paddle. This script loads the saved network and does not retrain.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
modelPath = fullfile(repoRoot, 'trained_networks', 'Pong_network.mat');
S = load(modelPath);
assert(isfield(S, 'split') && isfield(S.split, 'trainIdx') && isfield(S.split, 'valIdx'), 'Pong network lacks split metadata; rerun Pong.mlx.');
assert(isfield(S, 'metadata') && strcmp(S.metadata.normalization, 'none'), 'Pong network metadata must record no normalization.');
nObs = max([S.split.trainIdx(:); S.split.valIdx(:)]);
assertNoOverlap(S.split.trainIdx, S.split.valIdx, nObs, 'Pong');
net = loadSavedNetwork(modelPath);
rng(1, 'twister');

% Internal truncated networks used only for optional activation traces.
net2 = dlnetwork(net.Layers(1:5));
net3 = dlnetwork(net.Layers(1:3));
idx_neuron = [14068, 3042, 8509, 15830, 921]; %#ok<NASGU>

ballSpeed = 0.02;
maxBounceAngle = pi/4;
dt = 0.01;

paddleHeight = 0.2;
paddleWidth = 0.02;
ballRadius = 0.02;
oppPaddleX = 0.05;
netPaddleX = 1 - oppPaddleX - paddleWidth;

mirrorScore = 0;
netContactCount = 0;
netMissCount = 0;
maxNetworkHits = 100;
maxRounds = 500;
round = 0;

fig = figure('NumberTitle', 'off');
ax = axes('Parent', fig);
axis(ax, [0 1 0 1]);
axis(ax, 'equal');
set(ax, 'XTick', [], 'YTick', []);
box(ax, 'on');
hold(ax, 'on');

oppPaddleRect = rectangle(ax, 'Position', [oppPaddleX, 0.5 - paddleHeight/2, paddleWidth, paddleHeight], 'FaceColor', 'b');
netPaddleRect = rectangle(ax, 'Position', [netPaddleX, 0.5 - paddleHeight/2, paddleWidth, paddleHeight], 'FaceColor', 'g');
ballCircle = rectangle(ax, 'Position', [0.5-ballRadius, 0.5-ballRadius, 2*ballRadius, 2*ballRadius], 'Curvature', [1, 1], 'FaceColor', 'r');
scoreText = text(ax, 0.5, 0.97, '', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

while netContactCount < maxNetworkHits && round < maxRounds && ishandle(fig)
    round = round + 1;
    ballPos = [0.5, 0.5];
    theta = (pi/4)*rand() + pi/8;
    if rand() < 0.5
        ballVel = ballSpeed * [cos(theta), sin(theta)];
    else
        ballVel = ballSpeed * [-cos(theta), sin(theta)];
    end

    oppPaddleY = 0.5;
    netPaddleY = 0.5;
    updateObjects(oppPaddleRect, netPaddleRect, ballCircle, oppPaddleX, netPaddleX, oppPaddleY, netPaddleY, ballPos, paddleHeight, paddleWidth, ballRadius);
    updateScore(scoreText, round, mirrorScore, netContactCount, netMissCount);

    roundActive = true;
    cnt1 = 1;
    outy = [];
    truey = [];
    A2 = [];
    A1 = [];

    while roundActive && ishandle(fig)
        oppPaddleY = min(max(ballPos(2), paddleHeight/2), 1 - paddleHeight/2);

        state = [ballPos, ballVel, netPaddleY];
        netAction = scoresToAction(predict(net, state));
        if netAction == "Up"
            netPaddleY = max(netPaddleY - 0.025, paddleHeight/2);
        elseif netAction == "Down"
            netPaddleY = min(netPaddleY + 0.025, 1 - paddleHeight/2);
        end

        ballPos = ballPos + ballVel;
        if ballPos(2) - ballRadius <= 0 || ballPos(2) + ballRadius >= 1
            ballVel(2) = -ballVel(2);
        end

        A2(cnt1, :) = numericOutput(predict(net2, state))(:).'; %#ok<SAGROW>
        A1(cnt1, :) = numericOutput(predict(net3, state))(:).'; %#ok<SAGROW>
        outy(cnt1) = netPaddleY; %#ok<SAGROW>
        truey(cnt1, :) = ballPos; %#ok<SAGROW>
        cnt1 = cnt1 + 1;

        if ballPos(1) - ballRadius <= oppPaddleX + paddleWidth
            if ballPos(2) >= oppPaddleY - paddleHeight/1.5 && ballPos(2) <= oppPaddleY + paddleHeight/1.5
                relIntersect = (ballPos(2) - oppPaddleY) / (paddleHeight/2);
                angle = relIntersect * maxBounceAngle;
                ballVel(1) = abs(ballSpeed * cos(angle));
                ballVel(2) = ballSpeed * sin(angle);
                ballPos(1) = oppPaddleX + paddleWidth + ballRadius;
            else
                mirrorScore = mirrorScore + 1;
                roundActive = false;
            end
        end

        if ballPos(1) + ballRadius >= netPaddleX
            if ballPos(2) >= netPaddleY - paddleHeight/1.5 && ballPos(2) <= netPaddleY + paddleHeight/1.5
                netContactCount = netContactCount + 1;
                relIntersect = (ballPos(2) - netPaddleY) / (paddleHeight/2);
                angle = relIntersect * maxBounceAngle;
                ballVel(1) = -abs(ballSpeed * cos(angle));
                ballVel(2) = ballSpeed * sin(angle);
                ballPos(1) = netPaddleX - ballRadius;
            else
                netMissCount = netMissCount + 1;
                roundActive = false;
            end
        end

        if ballPos(1) < 0 || ballPos(1) > 1
            roundActive = false;
        end

        updateObjects(oppPaddleRect, netPaddleRect, ballCircle, oppPaddleX, netPaddleX, oppPaddleY, netPaddleY, ballPos, paddleHeight, paddleWidth, ballRadius);
        updateScore(scoreText, round, mirrorScore, netContactCount, netMissCount);
        drawnow;
        pause(dt);
    end
end

assert(netContactCount > 0 || netMissCount > 0, 'Pong simulation did not reach the network paddle.');
fprintf('After %d rounds:\\n', round);
fprintf('Mirror (Left) Score: %d\\n', mirrorScore);
fprintf('Network Hits: %d\\n', netContactCount);
fprintf('Network Misses: %d\\n', netMissCount);

${common}

function action = scoresToAction(scores)
    scores = numericOutput(scores);
    scores = scores(:);
    [~, idx] = max(scores);
    actions = ["Up", "Down", "Stay"];
    action = actions(idx);
end

function updateObjects(oppPaddleRect, netPaddleRect, ballCircle, oppPaddleX, netPaddleX, oppPaddleY, netPaddleY, ballPos, paddleHeight, paddleWidth, ballRadius)
    set(oppPaddleRect, 'Position', [oppPaddleX, oppPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
    set(netPaddleRect, 'Position', [netPaddleX, netPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
    set(ballCircle, 'Position', [ballPos(1)-ballRadius, ballPos(2)-ballRadius, 2*ballRadius, 2*ballRadius]);
end

function updateScore(scoreText, round, mirrorScore, netContactCount, netMissCount)
    scoreText.String = sprintf('Round %d | Mirror: %d | Network hits: %d | Network misses: %d', round, mirrorScore, netContactCount, netMissCount);
end

function assertNoOverlap(trainIdx, valIdx, n, label)
    assert(~isempty(trainIdx) && ~isempty(valIdx), '%s train/validation split is empty.', label);
    assert(all(trainIdx >= 1 & trainIdx <= n) && all(valIdx >= 1 & valIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(numel(unique([trainIdx(:); valIdx(:)])) == n, '%s duplicate or missing split indices.', label);
end
`;

makeMlx(path.join(repo, "examples", "control", "motor_control.mlx"), motorControlTrain);
makeMlx(path.join(repo, "examples", "pong", "Pong.mlx"), pongTrain);
makeMlx(path.join(repo, "tests", "test_pong_saved_network.mlx"), pongTest);
makeMlx(path.join(repo, "tests", "Pong_test.mlx"), pongSimulationTest);

const lqrTest = `%% Test saved LQR network with hard leakage checks
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
modelPath = fullfile(repoRoot, 'trained_networks', 'LQR_network.mat');
S = load(modelPath);
requiredFields = {'net', 'mu_X', 'sigma_X', 'mu_Y', 'sigma_Y', 'split', 'testSequence'};
for i = 1:numel(requiredFields)
    assert(isfield(S, requiredFields{i}), 'LQR network lacks %s; rerun motor_control.mlx.', requiredFields{i});
end
assertEpisodeSplit(S.split.trainEps, S.split.valEps, 250, 'LQR');
assert(isfield(S.testSequence, 'trainTargets') && isfield(S.testSequence, 'testTargets'), 'LQR testSequence lacks train/test target metadata.');
assertNoSharedTargets(S.testSequence.trainTargets, S.testSequence.testTargets, 'LQR train/test target leakage detected.');
net = loadSavedNetwork(modelPath);
assertNetworkIO(net, numel(S.mu_X), 2, 'LQR');
fprintf('LQR saved network split, normalization, and target leakage checks passed.\\n');

${common}
${lqrHelpers()}

function assertNetworkIO(net, expectedInputs, expectedOutputs, label)
    y = numericOutput(predict(net, dlarray(zeros(expectedInputs, 1, 'single'), 'CB')));
    assert(numel(y) == expectedOutputs, '%s output size mismatch.', label);
end
`;

makeMlx(path.join(repo, "tests", "test_lqr_saved_network.mlx"), lqrTest);

console.log("Created leakage-safe legacy train/test live scripts.");
