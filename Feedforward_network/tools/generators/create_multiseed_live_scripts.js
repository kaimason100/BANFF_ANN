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
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-seed-"));
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
  childProcess.execFileSync("zip", ["-qr", outPath, "."], { cwd: work });
  fs.rmSync(work, { recursive: true, force: true });
}

const common = `
function repoRoot = locateRepoRoot()
% Return the release root when this live script is run from any folder.
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
end

function setWeightSeed(seedValue)
% Set the global offset used by customWeights*_git. Default scripts use zero.
    setappdata(0, 'NAR_WEIGHT_SEED_BASE', seedValue);
end

function clearWeightSeed()
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        rmappdata(0, 'NAR_WEIGHT_SEED_BASE');
    end
end

function assertSeedMetadata(S, task, seedValue, seedList)
    assert(isfield(S, 'metadata'), '%s seed %d lacks metadata.', task, seedValue);
    assert(isfield(S.metadata, 'task') && strcmp(char(S.metadata.task), char(task)), ...
        '%s seed %d metadata task does not match the tested task.', task, seedValue);
    assert(isfield(S.metadata, 'seed') && isequal(double(S.metadata.seed), double(seedValue)), ...
        '%s seed %d metadata seed does not match the file label.', task, seedValue);
    assert(isfield(S.metadata, 'seedList') && isequal(double(S.metadata.seedList(:)), double(seedList(:))), ...
        '%s seed %d metadata seed list does not match this test script.', task, seedValue);
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

function assertFeatureStatsFromTraining(stats, X, trainIdx, label)
    expected = fitFeatureStandardization(X(:, trainIdx));
    assertStatsClose(stats.mu, expected.mu, label, 'feature mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'feature standard deviation');
    assertStatsClose(stats.scale, expected.scale, label, 'feature scale');
end

function assertImageStatsFromTraining(stats, X, trainIdx, label)
    expected = fitImageStandardization(X(:, :, :, trainIdx));
    assertStatsClose(stats.mu, expected.mu, label, 'image mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'image standard deviation');
    assertStatsClose(stats.scale, expected.scale, label, 'image scale');
end

function assertStatsClose(actual, expected, label, statName)
    actual = double(actual);
    expected = double(expected);
    assert(isequal(size(actual), size(expected)), '%s %s size mismatch.', label, statName);
    tol = 1e-10 * max(1, max(abs(expected(:))));
    assert(all(abs(actual(:) - expected(:)) <= tol), '%s %s was not fitted from the training split.', label, statName);
end

function net = loadSavedNetwork(modelPath)
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
`;

const classificationTrain = `%% Train classification tasks across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');

tasks = ["iris", "breast_cancer", "car_quality", "mushroom"];
for task = tasks
    task = char(task);
    [X, labels, batchSize] = classificationDataset(repoRoot, task);
    [trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), 42, task);
    featureStats = fitFeatureStandardization(X(:, trainIdx));
    assertFeatureStatsFromTraining(featureStats, X, trainIdx, task);
    for seedValue = SEEDS
        fprintf('\\n[%s] Training seed %d\\n', task, seedValue);
        setWeightSeed(seedValue);
        try
        XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = labels(trainIdx);
        XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = labels(valIdx);
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
            'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', 0);
        [net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrain, net, 'crossentropy', options);
        split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
        metadata.task = task; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
        outDir = fullfile(OUT_ROOT, task); if ~exist(outDir, 'dir'), mkdir(outDir); end
        save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', task, seedValue)), 'net', 'info', 'split', 'featureStats', 'metadata', '-v7.3');
        clearWeightSeed();
        catch ME
            clearWeightSeed();
            rethrow(ME);
        end
    end
end

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

function [X, labels, batchSize] = classificationDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'iris'
            load fisheriris
            labels = categorical(species);
            X = meas.';
            batchSize = 150;
        case 'breast_cancer'
            S = load(fullfile(dataDir, 'breast_cancer_dataset.mat')); [labels, X] = breastCancerData(S.data); batchSize = 256;
        case 'car_quality'
            S = load(fullfile(dataDir, 'car_dataset.mat')); D = S.data; labels = categorical(D(:, 7)); D(~isfinite(D)) = 0; X = D(:, 1:6).'; batchSize = 256;
        case 'mushroom'
            S = load(fullfile(dataDir, 'mushroom_dataset.mat')); D = S.data; labels = categorical(D(:, 1)); D(~isfinite(D)) = 0; X = D(:, 2:end).'; batchSize = 4096;
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
`;

const classificationTest = `%% Aggregate five-seed classification test results
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
tasks = ["iris", "breast_cancer", "car_quality", "mushroom"];
for task = tasks
    task = char(task);
    [X, labels, ~] = classificationDataset(repoRoot, task);
    [expectedTrainIdx, expectedValIdx, expectedTestIdx] = fixedSplit(size(X, 2), 42, task);
    acc = nan(numel(SEEDS), 1);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', task, sprintf('%s_seed_%03d_network.mat', task, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assert(isfield(S, 'split'), 'Seeded model lacks split metadata: %s', modelPath);
        assertSeedMetadata(S, task, seedValue, SEEDS);
        assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
        assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), task);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
            '%s seed %d split metadata does not match the deterministic held-out test split.', task, seedValue);
        assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, task);
        net = loadSavedNetwork(modelPath);
        XTest = applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats); YTest = labels(S.split.testIdx);
        YTrain = labels(S.split.trainIdx);
        scores = numericOutput(predict(net, dlarray(single(XTest), 'CB')));
        [~, idx] = max(scores, [], 1);
        yHat = categorical(idx(:), 1:numel(categories(YTrain)), categories(YTrain));
        acc(i) = mean(yHat(:) == YTest(:));
        fprintf('%s seed %d accuracy: %.2f%%\\n', task, seedValue, 100*acc(i));
    end
    fprintf('%s five-seed accuracy mean ± SD: %.2f ± %.2f%%\\n', task, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));
end

${common}
${classificationTrain.slice(classificationTrain.indexOf('function [trainIdx'))}
`;

const regressionTrain = `%% Train regression tasks across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');

tasks = ["abalone", "toyota"];
for task = tasks
    task = char(task);
    [X, y, batchSize] = regressionDataset(repoRoot, task);
    [trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), task);
    featureStats = fitFeatureStandardization(X(:, trainIdx));
    assertFeatureStatsFromTraining(featureStats, X, trainIdx, task);
    yMu = mean(y(trainIdx)); yStd = std(y(trainIdx)); if yStd == 0, yStd = 1; end
    for seedValue = SEEDS
        fprintf('\\n[%s] Training seed %d\\n', task, seedValue);
        setWeightSeed(seedValue);
        try
        XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = (y(trainIdx) - yMu) / yStd;
        XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = (y(valIdx) - yMu) / yStd;
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
            'ValidationData', {dlarray(single(XVal), 'CB'), YVal}, ...
            'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', 0);
        [net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrain, net, 'mse', options);
        split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
        targetStats.yMu = yMu; targetStats.yStd = yStd;
        metadata.task = task; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
        outDir = fullfile(OUT_ROOT, task); if ~exist(outDir, 'dir'), mkdir(outDir); end
        save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', task, seedValue)), 'net', 'info', 'split', 'featureStats', 'targetStats', 'metadata', '-v7.3');
        clearWeightSeed();
        catch ME
            clearWeightSeed();
            rethrow(ME);
        end
    end
end

${common}

function [trainIdx, valIdx, testIdx] = fixedSplit(n, task)
    if strcmp(task, 'abalone'), rng(1, 'twister'); else, rng(42, 'twister'); end
    cv = cvpartition(n, 'HoldOut', 0.2);
    trainFullIdx = find(training(cv));
    testIdx = find(test(cv));
    cvVal = cvpartition(numel(trainFullIdx), 'HoldOut', 0.2);
    trainIdx = trainFullIdx(training(cvVal));
    valIdx = trainFullIdx(test(cvVal));
    assertDisjointSplits(trainIdx, valIdx, testIdx, n, task);
end

function [X, y, batchSize] = regressionDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'abalone'
            S = load(fullfile(dataDir, 'abalone_dataset.mat')); D = S.data; D(~isfinite(D)) = 0; y = D(:, end); X = D(:, 1:end-1).'; batchSize = size(X, 2);
        case 'toyota'
            S = load(fullfile(dataDir, 'toyota_dataset.mat')); D = S.data; D(~isfinite(D)) = 0; y = D(:, 3); X = D(:, [1, 2, 4:size(D, 2)]).'; batchSize = min(4096, size(X, 2));
        otherwise
            error('Unknown regression task: %s', task);
    end
    X(~isfinite(X)) = 0;
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
`;

const regressionTest = `%% Aggregate five-seed regression test results
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
tasks = ["abalone", "toyota"];
for task = tasks
    task = char(task);
    [X, y, ~] = regressionDataset(repoRoot, task);
    [expectedTrainIdx, expectedValIdx, expectedTestIdx] = fixedSplit(size(X, 2), task);
    metrics = struct('meanError', [], 'sdError', [], 'rmse', [], 'r', [], 'p', [], 'medAE', []);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', task, sprintf('%s_seed_%03d_network.mat', task, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assertSeedMetadata(S, task, seedValue, SEEDS);
        assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
        assert(isfield(S, 'targetStats') && isfield(S.targetStats, 'yMu') && isfield(S.targetStats, 'yStd'), 'Seeded regression model lacks target normalization metadata: %s', modelPath);
        assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), task);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
            '%s seed %d split metadata does not match the deterministic held-out test split.', task, seedValue);
        assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, task);
        net = loadSavedNetwork(modelPath);
        XTest = applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats);
        YpredN = numericOutput(predict(net, dlarray(single(XTest), 'CB')));
        Ypred = YpredN(:) * S.targetStats.yStd + S.targetStats.yMu;
        Ytrue = double(y(S.split.testIdx));
        err = Ypred(:) - Ytrue(:);
        metrics.meanError(i,1) = mean(err, 'omitnan');
        metrics.sdError(i,1) = std(err, 0, 'omitnan');
        metrics.rmse(i,1) = sqrt(mean(err.^2, 'omitnan'));
        metrics.medAE(i,1) = median(abs(err), 'omitnan');
        [metrics.r(i,1), metrics.p(i,1)] = corrStats(Ypred(:), Ytrue(:));
        fprintf('%s seed %d: meanErr=%.4f, sdErr=%.4f, RMSE=%.4f, r=%.4f, p=%.4g\\n', task, seedValue, metrics.meanError(i), metrics.sdError(i), metrics.rmse(i), metrics.r(i), metrics.p(i));
    end
    fprintf('%s five-seed RMSE mean ± SD: %.4f ± %.4f\\n', task, mean(metrics.rmse, 'omitnan'), std(metrics.rmse, 0, 'omitnan'));
    fprintf('%s five-seed signed error mean ± SD across seeds: %.4f ± %.4f\\n', task, mean(metrics.meanError, 'omitnan'), std(metrics.meanError, 0, 'omitnan'));
    figure('Color', 'w'); histogram(metrics.rmse); grid on; xlabel('RMSE'); ylabel('Seed count'); title(sprintf('%s five-seed RMSE distribution', task));
end

${common}
${regressionTrain.slice(regressionTrain.indexOf('function [trainIdx'))}

function [r, p] = corrStats(yPred, yTrue)
    valid = isfinite(yPred) & isfinite(yTrue);
    yPred = yPred(valid); yTrue = yTrue(valid);
    if numel(yPred) < 3 || std(yPred) == 0 || std(yTrue) == 0
        r = NaN; p = NaN; return
    end
    [R, P] = corrcoef(yPred, yTrue); r = R(1, 2); p = P(1, 2);
end
`;

const classificationHelpers = classificationTrain.slice(classificationTrain.indexOf("function [trainIdx"));
const regressionHelpers = regressionTrain.slice(regressionTrain.indexOf("function [trainIdx"));

function classificationSingleTrain(task) {
  return `%% Train ${task} classification across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');
TASK = '${task}';

[X, labels, batchSize] = classificationDataset(repoRoot, TASK);
[trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), 42, TASK);
featureStats = fitFeatureStandardization(X(:, trainIdx));
assertFeatureStatsFromTraining(featureStats, X, trainIdx, TASK);
for seedValue = SEEDS
    fprintf('\\n[%s] Training seed %d\\n', TASK, seedValue);
    setWeightSeed(seedValue);
    try
    XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = labels(trainIdx);
    XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = labels(valIdx);
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
        'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', 0);
    [net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrain, net, 'crossentropy', options);
    split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
    metadata.task = TASK; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
    outDir = fullfile(OUT_ROOT, TASK); if ~exist(outDir, 'dir'), mkdir(outDir); end
    save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', TASK, seedValue)), 'net', 'info', 'split', 'featureStats', 'metadata', '-v7.3');
    clearWeightSeed();
    catch ME
        clearWeightSeed();
        rethrow(ME);
    end
end

${common}
${classificationHelpers}
`;
}

function classificationSingleTest(task) {
  return `%% Test ${task} classification across five saved weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
TASK = '${task}';

[X, labels, ~] = classificationDataset(repoRoot, TASK);
[expectedTrainIdx, expectedValIdx, expectedTestIdx] = fixedSplit(size(X, 2), 42, TASK);
acc = nan(numel(SEEDS), 1);
for i = 1:numel(SEEDS)
    seedValue = SEEDS(i);
    modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', TASK, sprintf('%s_seed_%03d_network.mat', TASK, seedValue));
    assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
    S = load(modelPath);
    assert(isfield(S, 'split'), 'Seeded model lacks split metadata: %s', modelPath);
    assertSeedMetadata(S, TASK, seedValue, SEEDS);
    assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
    assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), TASK);
    assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
        '%s seed %d split metadata does not match the deterministic held-out test split.', TASK, seedValue);
    assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, TASK);
    net = loadSavedNetwork(modelPath);
    XTest = applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats); YTest = labels(S.split.testIdx);
    YTrain = labels(S.split.trainIdx);
    scores = numericOutput(predict(net, dlarray(single(XTest), 'CB')));
    [~, idx] = max(scores, [], 1);
    yHat = categorical(idx(:), 1:numel(categories(YTrain)), categories(YTrain));
    acc(i) = mean(yHat(:) == YTest(:));
    fprintf('%s seed %d accuracy: %.2f%%\\n', TASK, seedValue, 100*acc(i));
end
fprintf('%s five-seed accuracy mean ± SD: %.2f ± %.2f%%\\n', TASK, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));

${common}
${classificationHelpers}
`;
}

function regressionSingleTrain(task) {
  return `%% Train ${task} regression across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
MAX_EPOCHS = 10000;
LR = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');
TASK = '${task}';

[X, y, batchSize] = regressionDataset(repoRoot, TASK);
[trainIdx, valIdx, testIdx] = fixedSplit(size(X, 2), TASK);
featureStats = fitFeatureStandardization(X(:, trainIdx));
assertFeatureStatsFromTraining(featureStats, X, trainIdx, TASK);
yMu = mean(y(trainIdx)); yStd = std(y(trainIdx)); if yStd == 0, yStd = 1; end
for seedValue = SEEDS
    fprintf('\\n[%s] Training seed %d\\n', TASK, seedValue);
    setWeightSeed(seedValue);
    try
    XTrain = applyFeatureStandardization(X(:, trainIdx), featureStats); YTrain = (y(trainIdx) - yMu) / yStd;
    XVal = applyFeatureStandardization(X(:, valIdx), featureStats); YVal = (y(valIdx) - yMu) / yStd;
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
        'ValidationData', {dlarray(single(XVal), 'CB'), YVal}, ...
        'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', 0);
    [net, info] = trainnet(dlarray(single(XTrain), 'CB'), YTrain, net, 'mse', options);
    split.trainIdx = trainIdx; split.valIdx = valIdx; split.testIdx = testIdx;
    targetStats.yMu = yMu; targetStats.yStd = yStd;
    metadata.task = TASK; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
    outDir = fullfile(OUT_ROOT, TASK); if ~exist(outDir, 'dir'), mkdir(outDir); end
    save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', TASK, seedValue)), 'net', 'info', 'split', 'featureStats', 'targetStats', 'metadata', '-v7.3');
    clearWeightSeed();
    catch ME
        clearWeightSeed();
        rethrow(ME);
    end
end

${common}
${regressionHelpers}
`;
}

function regressionSingleTest(task) {
  return `%% Test ${task} regression across five saved weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
TASK = '${task}';

[X, y, ~] = regressionDataset(repoRoot, TASK);
[expectedTrainIdx, expectedValIdx, expectedTestIdx] = fixedSplit(size(X, 2), TASK);
metrics = struct('meanError', [], 'sdError', [], 'rmse', [], 'r', [], 'p', [], 'medAE', []);
for i = 1:numel(SEEDS)
    seedValue = SEEDS(i);
    modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', TASK, sprintf('%s_seed_%03d_network.mat', TASK, seedValue));
    assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
    S = load(modelPath);
    assert(isfield(S, 'split'), 'Seeded model lacks split metadata: %s', modelPath);
    assertSeedMetadata(S, TASK, seedValue, SEEDS);
    assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
    assert(isfield(S, 'targetStats') && isfield(S.targetStats, 'yMu') && isfield(S.targetStats, 'yStd'), 'Seeded regression model lacks target normalization metadata: %s', modelPath);
    assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), TASK);
    assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
        '%s seed %d split metadata does not match the deterministic held-out test split.', TASK, seedValue);
    assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, TASK);
    net = loadSavedNetwork(modelPath);
    XTest = applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats);
    YpredN = numericOutput(predict(net, dlarray(single(XTest), 'CB')));
    Ypred = YpredN(:) * S.targetStats.yStd + S.targetStats.yMu;
    Ytrue = double(y(S.split.testIdx));
    err = Ypred(:) - Ytrue(:);
    metrics.meanError(i,1) = mean(err, 'omitnan');
    metrics.sdError(i,1) = std(err, 0, 'omitnan');
    metrics.rmse(i,1) = sqrt(mean(err.^2, 'omitnan'));
    metrics.medAE(i,1) = median(abs(err), 'omitnan');
    [metrics.r(i,1), metrics.p(i,1)] = corrStats(Ypred(:), Ytrue(:));
    fprintf('%s seed %d: meanErr=%.4f, sdErr=%.4f, RMSE=%.4f, r=%.4f, p=%.4g\\n', TASK, seedValue, metrics.meanError(i), metrics.sdError(i), metrics.rmse(i), metrics.r(i), metrics.p(i));
end
fprintf('%s five-seed RMSE mean ± SD: %.4f ± %.4f\\n', TASK, mean(metrics.rmse, 'omitnan'), std(metrics.rmse, 0, 'omitnan'));
fprintf('%s five-seed signed error mean ± SD across seeds: %.4f ± %.4f\\n', TASK, mean(metrics.meanError, 'omitnan'), std(metrics.meanError, 0, 'omitnan'));
figure('Color', 'w'); histogram(metrics.rmse); grid on; xlabel('RMSE'); ylabel('Seed count'); title(sprintf('%s five-seed RMSE distribution', TASK));

${common}
${regressionHelpers}

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

makeMlx(path.join(repo, "examples", "classification", "train_classification_5seeds.mlx"), classificationTrain);
makeMlx(path.join(repo, "tests", "test_classification_5seeds.mlx"), classificationTest);
makeMlx(path.join(repo, "examples", "regression", "train_regression_5seeds.mlx"), regressionTrain);
makeMlx(path.join(repo, "tests", "test_regression_5seeds.mlx"), regressionTest);

for (const task of ["iris", "breast_cancer", "car_quality", "mushroom"]) {
  makeMlx(path.join(repo, "examples", "classification", `train_${task}_5seeds.mlx`), classificationSingleTrain(task));
  makeMlx(path.join(repo, "tests", `test_${task}_5seeds.mlx`), classificationSingleTest(task));
}

for (const task of ["abalone", "toyota"]) {
  makeMlx(path.join(repo, "examples", "regression", `train_${task}_5seeds.mlx`), regressionSingleTrain(task));
  makeMlx(path.join(repo, "tests", `test_${task}_5seeds.mlx`), regressionSingleTest(task));
}

const mnistTrain = `%% Train MNIST-family classification tasks across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
NUM_EPOCHS = 1000;
BATCH_SIZE = 500;
LEARN_RATE = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');
DATA_DIR = fullfile(repoRoot, 'data');

targetList = ["afro_mnist_ethiopic.mat", "afro_mnist_nko.mat", "afro_mnist_osmanya.mat", ...
    "afro_mnist_vai.mat", "kmnist.mat", "mnist_fashion.mat", "mnist.mat"];

for ds = targetList
    ds = char(ds);
    dataPath = fullfile(DATA_DIR, ds);
    assert(isfile(dataPath), 'Missing MNIST-family dataset: %s', dataPath);
    D = load(dataPath);
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, ds);
    [trainIdx, valIdx] = mnistTrainValSplit(numel(YTrainFull), 42, ds);
    imageStats = fitImageStandardization(XTrainFull(:, :, :, trainIdx));
    assertImageStatsFromTraining(imageStats, XTrainFull, trainIdx, ds);
    outName = erase(ds, '.mat');
    for seedValue = SEEDS
        fprintf('\\n[%s] Training seed %d\\n', outName, seedValue);
        setWeightSeed(seedValue);
        try
        XTrain = applyImageStandardization(XTrainFull(:, :, :, trainIdx), imageStats); YTrain = YTrainFull(trainIdx);
        XVal = applyImageStandardization(XTrainFull(:, :, :, valIdx), imageStats); YVal = YTrainFull(valIdx);
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
        XTestEval = applyImageStandardization(XTest, imageStats);
        scores = numericOutput(predictImageBatches(net, XTestEval, 256));
        predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
        yHat = categorical(predIdx(:), 1:numOut);
        acc = 100 * mean(yHat(:) == YTest(:));
        split.trainIdx = trainIdx; split.valIdx = valIdx; split.testSplit = 'official_test_struct';
        metadata.task = outName; metadata.dataset = ds; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
        outDir = fullfile(OUT_ROOT, outName); if ~exist(outDir, 'dir'), mkdir(outDir); end
        save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', outName, seedValue)), 'net', 'info', 'acc', 'split', 'imageStats', 'metadata', '-v7.3');
        clearWeightSeed();
        fprintf('[%s] seed %d test accuracy: %.2f%%\\n', outName, seedValue, acc);
        catch ME
            clearWeightSeed();
            rethrow(ME);
        end
    end
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

const mnistTest = `%% Aggregate five-seed MNIST-family classification test results
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
DATA_DIR = fullfile(repoRoot, 'data');
targetList = ["afro_mnist_ethiopic.mat", "afro_mnist_nko.mat", "afro_mnist_osmanya.mat", ...
    "afro_mnist_vai.mat", "kmnist.mat", "mnist_fashion.mat", "mnist.mat"];

for ds = targetList
    ds = char(ds);
    outName = erase(ds, '.mat');
    D = load(fullfile(DATA_DIR, ds));
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, ds);
    [expectedTrainIdx, expectedValIdx] = mnistTrainValSplit(numel(YTrainFull), 42, ds);
    acc = nan(numel(SEEDS), 1);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', outName, sprintf('%s_seed_%03d_network.mat', outName, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assertSeedMetadata(S, outName, seedValue, SEEDS);
        assert(isfield(S, 'imageStats'), 'MNIST seeded model lacks train-fitted image normalization metadata.');
        assert(isfield(S, 'split') && strcmp(S.split.testSplit, 'official_test_struct'), 'MNIST seeded model lacks test-split metadata.');
        assert(isfield(S.split, 'trainIdx') && isfield(S.split, 'valIdx'), 'MNIST seeded model lacks train/validation split metadata.');
        assert(isempty(intersect(S.split.trainIdx(:), S.split.valIdx(:))), '%s seed %d train/validation leakage.', outName, seedValue);
        assert(all(S.split.trainIdx(:) >= 1 & S.split.trainIdx(:) <= numel(YTrainFull)) && all(S.split.valIdx(:) >= 1 & S.split.valIdx(:) <= numel(YTrainFull)), ...
            '%s seed %d train/validation split indices are out of range.', outName, seedValue);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)), ...
            '%s seed %d train/validation split metadata does not match the deterministic training split.', outName, seedValue);
        assertImageStatsFromTraining(S.imageStats, XTrainFull, S.split.trainIdx, outName);
        net = loadSavedNetwork(modelPath);
        XTestEval = applyImageStandardization(XTest, S.imageStats);
        scores = numericOutput(predictImageBatches(net, XTestEval, 256));
        predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
        yHat = categorical(predIdx(:), 1:numOut);
        acc(i) = mean(yHat(:) == YTest(:));
        fprintf('%s seed %d accuracy: %.2f%%\\n', outName, seedValue, 100*acc(i));
    end
    fprintf('%s five-seed accuracy mean ± SD: %.2f ± %.2f%%\\n', outName, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));
end

${common}
${mnistTrain.slice(mnistTrain.indexOf('function [XTrain'))}
`;

const mnistHelpers = mnistTrain.slice(mnistTrain.indexOf("function [XTrain"));

function mnistSingleTrain(dataset) {
  return `%% Train ${dataset.replace(".mat", "")} across five fixed weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
WIDTH = 16000;
NUM_EPOCHS = 1000;
BATCH_SIZE = 500;
LEARN_RATE = 0.01;
OUT_ROOT = fullfile(repoRoot, 'trained_networks', 'seeded');
DATA_DIR = fullfile(repoRoot, 'data');
DATASET = '${dataset}';

dataPath = fullfile(DATA_DIR, DATASET);
assert(isfile(dataPath), 'Missing MNIST-family dataset: %s', dataPath);
D = load(dataPath);
[XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, DATASET);
[trainIdx, valIdx] = mnistTrainValSplit(numel(YTrainFull), 42, DATASET);
imageStats = fitImageStandardization(XTrainFull(:, :, :, trainIdx));
assertImageStatsFromTraining(imageStats, XTrainFull, trainIdx, DATASET);
outName = erase(DATASET, '.mat');
for seedValue = SEEDS
    fprintf('\\n[%s] Training seed %d\\n', outName, seedValue);
    setWeightSeed(seedValue);
    try
    XTrain = applyImageStandardization(XTrainFull(:, :, :, trainIdx), imageStats); YTrain = YTrainFull(trainIdx);
    XVal = applyImageStandardization(XTrainFull(:, :, :, valIdx), imageStats); YVal = YTrainFull(valIdx);
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
    XTestEval = applyImageStandardization(XTest, imageStats);
    scores = numericOutput(predictImageBatches(net, XTestEval, 256));
    predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
    yHat = categorical(predIdx(:), 1:numOut);
    acc = 100 * mean(yHat(:) == YTest(:));
    split.trainIdx = trainIdx; split.valIdx = valIdx; split.testSplit = 'official_test_struct';
    metadata.task = outName; metadata.dataset = DATASET; metadata.seed = seedValue; metadata.seedList = SEEDS; metadata.width = WIDTH;
    outDir = fullfile(OUT_ROOT, outName); if ~exist(outDir, 'dir'), mkdir(outDir); end
    save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', outName, seedValue)), 'net', 'info', 'acc', 'split', 'imageStats', 'metadata', '-v7.3');
    clearWeightSeed();
    fprintf('[%s] seed %d test accuracy: %.2f%%\\n', outName, seedValue, acc);
    catch ME
        clearWeightSeed();
        rethrow(ME);
    end
end

${common}
${mnistHelpers}
`;
}

function mnistSingleTest(dataset) {
  return `%% Test ${dataset.replace(".mat", "")} across five saved weight seeds
repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

SEEDS = [101 202 303 404 505];
DATA_DIR = fullfile(repoRoot, 'data');
DATASET = '${dataset}';

outName = erase(DATASET, '.mat');
D = load(fullfile(DATA_DIR, DATASET));
[XTrainFull, YTrainFull, XTest, YTest, numOut] = mnistData(D.training, D.test, DATASET);
[expectedTrainIdx, expectedValIdx] = mnistTrainValSplit(numel(YTrainFull), 42, DATASET);
acc = nan(numel(SEEDS), 1);
for i = 1:numel(SEEDS)
    seedValue = SEEDS(i);
    modelPath = fullfile(repoRoot, 'trained_networks', 'seeded', outName, sprintf('%s_seed_%03d_network.mat', outName, seedValue));
    assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
    S = load(modelPath);
    assertSeedMetadata(S, outName, seedValue, SEEDS);
    assert(isfield(S, 'imageStats'), 'MNIST seeded model lacks train-fitted image normalization metadata.');
    assert(isfield(S, 'split') && strcmp(S.split.testSplit, 'official_test_struct'), 'MNIST seeded model lacks test-split metadata.');
    assert(isfield(S.split, 'trainIdx') && isfield(S.split, 'valIdx'), 'MNIST seeded model lacks train/validation split metadata.');
    assert(isempty(intersect(S.split.trainIdx(:), S.split.valIdx(:))), '%s seed %d train/validation leakage.', outName, seedValue);
    assert(all(S.split.trainIdx(:) >= 1 & S.split.trainIdx(:) <= numel(YTrainFull)) && all(S.split.valIdx(:) >= 1 & S.split.valIdx(:) <= numel(YTrainFull)), ...
        '%s seed %d train/validation split indices are out of range.', outName, seedValue);
    assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)), ...
        '%s seed %d train/validation split metadata does not match the deterministic training split.', outName, seedValue);
    assertImageStatsFromTraining(S.imageStats, XTrainFull, S.split.trainIdx, outName);
    net = loadSavedNetwork(modelPath);
    XTestEval = applyImageStandardization(XTest, S.imageStats);
    scores = numericOutput(predictImageBatches(net, XTestEval, 256));
    predIdx = scoresToClassIndex(scores, numOut, numel(YTest));
    yHat = categorical(predIdx(:), 1:numOut);
    acc(i) = mean(yHat(:) == YTest(:));
    fprintf('%s seed %d accuracy: %.2f%%\\n', outName, seedValue, 100*acc(i));
end
fprintf('%s five-seed accuracy mean ± SD: %.2f ± %.2f%%\\n', outName, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));

${common}
${mnistHelpers}
`;
}

makeMlx(path.join(repo, "examples", "classification", "train_mnist_family_5seeds.mlx"), mnistTrain);
makeMlx(path.join(repo, "tests", "test_mnist_family_5seeds.mlx"), mnistTest);

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
  makeMlx(path.join(repo, "examples", "classification", `train_${task}_5seeds.mlx`), mnistSingleTrain(dataset));
  makeMlx(path.join(repo, "tests", `test_${task}_5seeds.mlx`), mnistSingleTest(dataset));
}

console.log("Created shared and per-task five-seed live scripts.");
