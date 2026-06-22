const fs = require("fs");
const path = require("path");
const os = require("os");
const childProcess = require("child_process");

const repo = process.cwd();
const template = path.join(repo, "examples", "classification", "Iris_classification.mlx");
const testsDir = path.join(repo, "tests");

function xmlEscapeCdata(text) {
  return text.replaceAll("]]>", "]]]]><![CDATA[>");
}

function documentXml(code) {
  return `<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:p><w:pPr><w:pStyle w:val="code"/></w:pPr><w:r><w:t><![CDATA[${xmlEscapeCdata(code)}]]></w:t></w:r></w:p></w:body></w:document>`;
}

function outputXml() {
  return `<?xml version="1.0" encoding="UTF-8"?><embeddedOutputs><metaData><evaluationState>manual</evaluationState><layoutState>code</layoutState><outputStatus>ready</outputStatus></metaData><outputArray type="array"/><regionArray type="array"/></embeddedOutputs>`;
}

function makeMlx(name, code) {
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-test-"));
  childProcess.execFileSync("unzip", ["-q", template, "-d", work]);
  fs.writeFileSync(path.join(work, "matlab", "document.xml"), documentXml(code));
  fs.writeFileSync(path.join(work, "matlab", "output.xml"), outputXml());
  const outPath = path.join(testsDir, name);
  fs.rmSync(outPath, { force: true });
  childProcess.execFileSync("zip", ["-qr", outPath, "."], { cwd: work });
  fs.rmSync(work, { recursive: true, force: true });
}

const sharedHelpers = `
function repoRoot = locateRepoRoot()
% Return the repository root when this live script is run from any folder.
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
            if isfolder(fullfile(candidate, 'trained_networks')) && ...
                    isfolder(fullfile(candidate, 'tests'))
                repoRoot = candidate;
                return
            end
            parent = fileparts(candidate);
            if strcmp(parent, candidate)
                break
            end
            candidate = parent;
        end
    end
    error('Could not locate repository root. Run from inside the release folder.');
end

function addReleasePaths(repoRoot)
% Add source and example folders needed by saved networks and helper code.
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(genpath(fullfile(repoRoot, 'examples')));
    addpath(genpath(fullfile(repoRoot, 'src', 'dynamical_systems')));
end

function net = loadSavedNetwork(modelPath)
% Load the first MATLAB network-like object from a saved training artifact.
    assert(isfile(modelPath), 'Missing saved network: %s', modelPath);
    S = load(modelPath);
    names = fieldnames(S);
    preferred = ["net", "dlnet", "trainedNet", "netObj"];
    for name = preferred
        fieldName = char(name);
        if isfield(S, fieldName)
            v = S.(fieldName);
            if isNetworkLike(v)
                net = v;
                return
            end
        end
    end
    for k = 1:numel(names)
        v = S.(names{k});
        if isNetworkLike(v)
            net = v;
            return
        end
    end
    error('No network-like variable found in %s', modelPath);
end

function tf = isNetworkLike(v)
% Accept common Deep Learning Toolbox network classes.
    c = lower(class(v));
    tf = contains(c, 'network') || strcmp(c, 'dlnetwork') || ...
         strcmp(c, 'seriesnetwork') || strcmp(c, 'dagnetwork');
end

function Y = predictCB(net, X)
% Predict with channel-by-batch data, falling back for older network classes.
    try
        Y = predict(net, dlarray(single(X), 'CB'));
    catch
        try
            Y = predict(net, single(X.'));
        catch
            Y = predict(net, single(X));
        end
    end
end

function A = numericOutput(Y)
% Convert dlarray, gpuArray, or ordinary numeric prediction output to double.
    if isa(Y, 'dlarray')
        A = extractdata(Y);
    else
        A = Y;
    end
    try
        A = gather(A);
    catch
    end
    A = double(A);
end

function pct = accuracyPercent(acc)
% Accept either fraction-style accuracy in [0,1] or already-percent accuracy.
    pct = double(acc);
    if isfinite(pct) && abs(pct) <= 1
        pct = 100 * pct;
    end
end

function assertDisjointSplits(trainIdx, valIdx, testIdx, n, label)
% Hard guard against train/validation/test index leakage.
    trainIdx = trainIdx(:);
    valIdx = valIdx(:);
    testIdx = testIdx(:);
    assert(~isempty(trainIdx), '%s train split is empty.', label);
    assert(~isempty(valIdx), '%s validation split is empty.', label);
    assert(~isempty(testIdx), '%s test split is empty.', label);
    allIdx = [trainIdx; valIdx; testIdx];
    assert(all(allIdx >= 1 & allIdx <= n), '%s split indices are out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage detected.', label);
    assert(isempty(intersect(trainIdx, testIdx)), '%s train/test leakage detected.', label);
    assert(isempty(intersect(valIdx, testIdx)), '%s validation/test leakage detected.', label);
    assert(numel(unique(allIdx)) == numel(allIdx), '%s duplicate split indices detected.', label);
end
`;

makeMlx("run_all_tests.mlx", `%% Run all saved-network release tests
% Execute each test live script. These tests load saved artifacts and do not
% retrain networks or reset task hyperparameters.

repoRoot = locateRepoRoot();
testDir = fullfile(repoRoot, 'tests');
tests = [
    "test_iris_saved_network.mlx"
    "test_tabular_saved_networks.mlx"
    "test_regression_saved_networks.mlx"
    "test_dynamical_system_saved_networks.mlx"
    "test_lqr_saved_network.mlx"
    "test_mnist_saved_networks.mlx"
    "test_pong_saved_network.mlx"
];

for i = 1:numel(tests)
    fprintf('\\n=== %s ===\\n', tests(i));
    run(fullfile(testDir, tests(i)));
end

fprintf('\\nAll available saved-network tests finished.\\n');

function repoRoot = locateRepoRoot()
% Return the repository root when this live script is run from any folder.
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
            if isfolder(fullfile(candidate, 'trained_networks')) && ...
                    isfolder(fullfile(candidate, 'tests'))
                repoRoot = candidate;
                return
            end
            parent = fileparts(candidate);
            if strcmp(parent, candidate)
                break
            end
            candidate = parent;
        end
    end
    error('Could not locate repository root. Run from inside the release folder.');
end
`);

makeMlx("test_iris_saved_network.mlx", `%% Test saved Iris network
% Loads trained_networks/Iris_network.mat and evaluates it on a reproducible Iris
% holdout split. No training is performed.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

net = loadSavedNetwork(fullfile(repoRoot, 'trained_networks', 'Iris_network.mat'));

load fisheriris
labels = categorical(species);
X = meas.';
X = (1/sqrt(4)) * ((X - mean(X, 2)) ./ std(X, 0, 2));

rng(42, 'twister');
cv = cvpartition(size(X, 2), 'HoldOut', 0.2);
trainFullIdx = find(training(cv));
testIdx = find(test(cv));
XTrainFull = X(:, trainFullIdx);
YTrainFull = labels(trainFullIdx);
XTest = X(:, testIdx);
YTest = labels(testIdx);

cvVal = cvpartition(sum(training(cv)), 'HoldOut', 0.2);
trainIdx = trainFullIdx(training(cvVal));
valIdx = trainFullIdx(test(cvVal));
assertDisjointSplits(trainIdx, valIdx, testIdx, size(X, 2), 'Iris');
YTrain = labels(trainIdx);

scores = numericOutput(predictCB(net, XTest));
[~, idx] = max(scores, [], 1);
predictedLabels = categorical(idx(:), 1:numel(categories(YTrain)), categories(YTrain));

acc = mean(predictedLabels(:) == YTest(:));
fprintf('Iris saved-network accuracy: %.2f%%\\n', 100*acc);

${sharedHelpers}`);

makeMlx("test_tabular_saved_networks.mlx", `%% Test saved tabular classification networks
% Evaluates saved classification networks using held-out test data only.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

tasks = struct( ...
    'label', {'Breast cancer', 'Car quality', 'Mushroom'}, ...
    'model', {'BC_network.mat', 'Car_quality_network.mat', 'Mushroom_network.mat'}, ...
    'dataset', {'breast_cancer_dataset.mat', 'car_dataset.mat', 'mushroom_dataset.mat'});

for t = 1:numel(tasks)
    task = tasks(t);
    net = loadSavedNetwork(fullfile(repoRoot, 'trained_networks', task.model));
    dataPath = fullfile(repoRoot, 'data', task.dataset);
    assert(isfile(dataPath), 'Missing required test dataset: %s', dataPath);

    S = load(dataPath);
    assert(isfield(S, 'data'), 'Dataset must contain variable data: %s', dataPath);
    [X, labels] = classificationData(task.dataset, S.data);

    rng(42, 'twister');
    cv = cvpartition(size(X, 2), 'HoldOut', 0.2);
    trainFullIdx = find(training(cv));
    testIdx = find(test(cv));
    XTrainFull = X(:, trainFullIdx);
    YTrainFull = labels(trainFullIdx);
    XTest = X(:, testIdx);
    YTest = labels(testIdx);

    cvVal = cvpartition(sum(training(cv)), 'HoldOut', 0.2);
    trainIdx = trainFullIdx(training(cvVal));
    valIdx = trainFullIdx(test(cvVal));
    assertDisjointSplits(trainIdx, valIdx, testIdx, size(X, 2), task.label);
    YTrain = labels(trainIdx);

    scores = numericOutput(predictCB(net, XTest));
    [~, idx] = max(scores, [], 1);
    predictedLabels = categorical(idx(:), 1:numel(categories(YTrain)), categories(YTrain));
    acc = mean(predictedLabels(:) == YTest(:));
    fprintf('%s saved-network accuracy: %.2f%%\\n', task.label, 100*acc);
end

${sharedHelpers}

function [X, labels] = classificationData(datasetName, data)
% Convert each prepared classification dataset to C-by-N inputs plus labels.
    switch datasetName
        case 'breast_cancer_dataset.mat'
            [labels, X] = extractBreastCancer(data);
        case 'car_dataset.mat'
            labels = categorical(data(:, 7));
            data(~isfinite(data)) = 0;
            X = data(:, 1:6).';
            C = size(X, 1);
            X = (1/sqrt(C)) * ((X - mean(X, 2)) ./ std(X, 0, 2));
        case 'mushroom_dataset.mat'
            labels = categorical(data(:, 1));
            data(~isfinite(data)) = 0;
            X = data(:, 2:end).';
            X = normalize(X, 2) / sqrt(size(X, 1));
        otherwise
            error('Unsupported classification dataset: %s', datasetName);
    end
    X(~isfinite(X)) = 0;
end

function [labels, X] = extractBreastCancer(data)
% Match the breast-cancer training script's supported data formats.
    if istable(data)
        labels = categorical(string(data{:, 2}));
        X = data{:, 3:end}.';
    elseif iscell(data)
        labels = categorical(string(data(:, 2)));
        X = cellfun(@double, data(:, 3:end)).';
    elseif isnumeric(data)
        y = data(:, 2);
        [~, ~, g] = unique(y(:));
        labels = categorical(g);
        X = data(:, 3:end).';
    elseif isstruct(data) && isfield(data, 'features') && isfield(data, 'labels')
        X = data.features.';
        labels = categorical(data.labels(:));
    else
        error('Unsupported breast-cancer data format.');
    end
    C = size(X, 1);
    X = (1/sqrt(C)) * ((X - mean(X, 2)) ./ std(X, 0, 2));
end
`);

makeMlx("test_regression_saved_networks.mlx", `%% Test saved regression networks
% Evaluates saved regression networks using held-out test data only.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

tasks = struct( ...
    'label', {'Abalone age', 'Toyota price'}, ...
    'model', {'Abalone_network.mat', 'Toyota_network.mat'}, ...
    'dataset', {'abalone_dataset.mat', 'toyota_dataset.mat'});

for t = 1:numel(tasks)
    task = tasks(t);
    net = loadSavedNetwork(fullfile(repoRoot, 'trained_networks', task.model));
    dataPath = fullfile(repoRoot, 'data', task.dataset);
    assert(isfile(dataPath), 'Missing required test dataset: %s', dataPath);

    S = load(dataPath);
    assert(isfield(S, 'data'), 'Dataset must contain variable data: %s', dataPath);
    [XTest, YTest, yMu, yStd] = regressionTestData(task.dataset, S.data);
    YpredNorm = numericOutput(predictCB(net, XTest));
    Ypred = YpredNorm(:) * yStd + yMu;
    Ytrue = double(YTest(:));
    signedError = Ypred(:) - Ytrue;
    absError = abs(signedError);
    rmse = sqrt(mean(signedError.^2, 'omitnan'));
    muErr = mean(signedError, 'omitnan');
    sdErr = std(signedError, 0, 'omitnan');
    [r, p] = correlationStats(Ypred(:), Ytrue);

    fprintf('%s signed error mean: %.4f\\n', task.label, muErr);
    fprintf('%s signed error SD: %.4f\\n', task.label, sdErr);
    fprintf('%s RMSE: %.4f\\n', task.label, rmse);
    fprintf('%s Pearson r: %.4f, p: %.4g\\n', task.label, r, p);
    fprintf('%s median absolute error: %.4f\\n', task.label, median(absError, 'omitnan'));

    figure('Color', 'w');
    histogram(signedError, 40);
    grid on;
    xlabel('Signed error');
    ylabel('Count');
    title(sprintf('%s signed error distribution', task.label));

    figure('Color', 'w');
    histogram(absError, 40);
    grid on;
    xlabel('Absolute error');
    ylabel('Count');
    title(sprintf('%s absolute error distribution', task.label));
end

${sharedHelpers}

function [XTest, YTest, yMu, yStd] = regressionTestData(datasetName, data)
% Rebuild the same held-out split and target scaling used by training.
    data(~isfinite(data)) = 0;
    switch datasetName
        case 'abalone_dataset.mat'
            y = data(:, end);
            X = data(:, 1:end-1).';
            X = (1/sqrt(size(X, 1))) * ((X - mean(X, 2)) ./ std(X, 0, 2));
            rng(1, 'twister');
        case 'toyota_dataset.mat'
            y = data(:, 3);
            X = data(:, [1, 2, 4:size(data, 2)]).';
            sig = std(X, 0, 2);
            sig(sig == 0) = 1;
            X = (1/sqrt(size(X, 1))) * ((X - mean(X, 2)) ./ sig);
            rng(42, 'twister');
        otherwise
            error('Unsupported regression dataset: %s', datasetName);
    end
    X(~isfinite(X)) = 0;
    cv = cvpartition(size(X, 2), 'HoldOut', 0.2);
    trainFullIdx = find(training(cv));
    testIdx = find(test(cv));
    XTrainFull = X(:, trainFullIdx);
    YTrainFull = y(trainFullIdx);
    XTest = X(:, testIdx);
    YTest = y(testIdx);
    cvVal = cvpartition(size(XTrainFull, 2), 'HoldOut', 0.2);
    trainIdx = trainFullIdx(training(cvVal));
    valIdx = trainFullIdx(test(cvVal));
    assertDisjointSplits(trainIdx, valIdx, testIdx, size(X, 2), datasetName);
    YTrain = y(trainIdx);
    yMu = mean(YTrain);
    yStd = std(YTrain);
    if yStd == 0
        yStd = 1;
    end
end

function [r, p] = correlationStats(yPred, yTrue)
% Return Pearson correlation and p-value after removing missing values.
    valid = isfinite(yPred) & isfinite(yTrue);
    yPred = yPred(valid);
    yTrue = yTrue(valid);
    if numel(yPred) < 3 || std(yPred) == 0 || std(yTrue) == 0
        r = NaN;
        p = NaN;
        return
    end
    [R, P] = corrcoef(yPred, yTrue);
    r = R(1, 2);
    p = P(1, 2);
end
`);


makeMlx("test_lqr_saved_network.mlx", `%% Test saved LQR two-link arm network
% Loads the saved motor-control network and evaluates it on fresh reaching
% targets using the same two-link arm dynamics and periodic error encoding as
% examples/control/motor_control.mlx. No training episodes are used as tests.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

networkPath = fullfile(repoRoot, 'trained_networks', 'LQR_network.mat');
net = loadSavedNetwork(networkPath);
Snet = load(networkPath);

[A, B, K, dt, T, t, l1, l2] = lqrArmSetup();
[mu_X, sigma_X, mu_Y, sigma_Y, numInputs, trainTargets] = lqrStatsOrRegenerate(Snet, A, B, K, dt, t, l1, l2);
assert(numInputs == 6, 'LQR network expects six encoded inputs.');
assertNetworkIO(net, numInputs, 2, 'LQR two-link arm');

numSims = 25;
rng(991, 'twister');
testTargets = generateLqrTargets(numSims, l1, l2);
assertNoSharedTargets(trainTargets, testTargets, 'LQR train/test target leakage detected.');

lqrFinishTimes = zeros(numSims, 1);
netFinishTimes = zeros(numSims, 1);
finalJointRmse = zeros(numSims, 1);
networkSuccessCount = 0;

for sim = 1:numSims
    xTargetJoint = inverseKinematicsTarget(testTargets(sim, :).', l1, l2);
    [xLqr, lqrFinishTimes(sim)] = simulateLqrArm(A, B, K, dt, t, xTargetJoint, testTargets(sim, :).', l1, l2);
    [xNet, netFinishTimes(sim), netFinished] = simulateNetworkArm(net, A, B, dt, t, xTargetJoint, testTargets(sim, :).', l1, l2, mu_X, sigma_X, mu_Y, sigma_Y, numInputs);
    networkSuccessCount = networkSuccessCount + double(netFinished);
    finalJointRmse(sim) = sqrt(mean((xNet(:, end) - xTargetJoint).^2));
    assert(all(isfinite(xLqr(:))) && all(isfinite(xNet(:))), 'Non-finite LQR motor-control trajectory.');
end

timeDiffs = netFinishTimes - lqrFinishTimes;
fprintf('LQR motor-control test simulations: %d\\n', numSims);
fprintf('Network successful simulations: %d\\n', networkSuccessCount);
fprintf('Mean finish-time difference (Network - LQR): %.4f s\\n', mean(timeDiffs, 'omitnan'));
fprintf('Final joint RMSE mean: %.4f\\n', mean(finalJointRmse, 'omitnan'));

assert(any(isfinite(finalJointRmse)), 'LQR motor-control test produced no finite RMSE values.');

${sharedHelpers}

function [A, B, K, dt, T, t, l1, l2] = lqrArmSetup()
% Shared two-link arm setup from the motor-control training script.
    A = [0 0 1 0;
         0 0 0 1;
         0 0 0 0;
         0 0 0 0];
    B = [0 0;
         0 0;
         1 0;
         0 1];
    Q = diag([100, 100, 1, 1]);
    R = diag([0.1, 0.1]);
    K = lqr(A, B, Q, R);
    dt = 0.01;
    T = 5;
    t = 0:dt:T;
    l1 = 1.0;
    l2 = 0.8;
end

function [mu_X, sigma_X, mu_Y, sigma_Y, numInputs, trainTargets] = lqrStatsOrRegenerate(Snet, A, B, K, dt, t, l1, l2)
% Prefer saved training normalisation stats; otherwise regenerate them exactly.
    if all(isfield(Snet, {'mu_X', 'sigma_X', 'mu_Y', 'sigma_Y'}))
        mu_X = Snet.mu_X;
        sigma_X = Snet.sigma_X;
        mu_Y = Snet.mu_Y;
        sigma_Y = Snet.sigma_Y;
        numInputs = numel(mu_X);
        [~, ~, trainTargets] = regenerateLqrTrainingStats(A, B, K, dt, t, l1, l2);
    else
        [stats, ~, trainTargets] = regenerateLqrTrainingStats(A, B, K, dt, t, l1, l2);
        mu_X = stats.mu_X;
        sigma_X = stats.sigma_X;
        mu_Y = stats.mu_Y;
        sigma_Y = stats.sigma_Y;
        numInputs = stats.numInputs;
    end
    sigma_X(sigma_X == 0) = 1;
    sigma_Y(sigma_Y == 0) = 1;
end

function [stats, splitInfo, trainTargets] = regenerateLqrTrainingStats(A, B, K, dt, t, l1, l2)
% Deterministically regenerate the supervised episodes used by motor_control.
    rng(1, 'twister');
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
            dist = endEffectorDistance(xSim(:, k+1), target, l1, l2);
            if dist < 0.01
                consecutiveCounter = consecutiveCounter + 1;
            else
                consecutiveCounter = 0;
            end
            if consecutiveCounter >= 10
                break
            end
        end
        inputsEpisodes{ep} = epiInputs;
        outputsEpisodes{ep} = epiOutputs;
    end
    idxEp = randperm(numEpisodes);
    numTrainEp = round(0.8 * numEpisodes);
    trainEps = idxEp(1:numTrainEp);
    valEps = idxEp(numTrainEp+1:end);
    assert(~isempty(trainEps) && ~isempty(valEps), 'LQR train/validation episode split is empty.');
    assert(isempty(intersect(trainEps, valEps)), 'LQR train/validation episode leakage detected.');
    assert(numel(unique([trainEps(:); valEps(:)])) == numEpisodes, 'LQR duplicate or missing episode split indices.');
    XTrainAll = localConcat(inputsEpisodes, trainEps);
    YTrainAll = localConcat(outputsEpisodes, trainEps);
    XValAll = localConcat(inputsEpisodes, valEps);
    assert(~isempty(XTrainAll) && ~isempty(XValAll), 'LQR regenerated train/validation samples are empty.');
    stats.mu_X = mean(XTrainAll, 2);
    stats.sigma_X = std(XTrainAll, 0, 2);
    stats.sigma_X(stats.sigma_X == 0) = 1;
    stats.mu_Y = mean(YTrainAll, 2);
    stats.sigma_Y = std(YTrainAll, 0, 2);
    stats.sigma_Y(stats.sigma_Y == 0) = 1;
    stats.numInputs = size(XTrainAll, 1);
    splitInfo.trainEps = trainEps;
    splitInfo.valEps = valEps;
    trainTargets = targetEE(trainEps, :);
end

function targets = generateLqrTargets(numTargets, l1, l2)
% Generate fresh reaching targets for held-out evaluation.
    targets = zeros(numTargets, 2);
    for i = 1:numTargets
        targets(i, :) = randomReachTarget(l1, l2).';
    end
end

function target = randomReachTarget(l1, l2)
% Random reachable end-effector target used by the original task.
    r = 0.5 + (l1 + l2 - 0.5) * rand;
    angle = 2*pi*rand;
    target = [r*cos(angle); r*sin(angle)];
end

function xTargetJoint = inverseKinematicsTarget(target, l1, l2)
% Convert end-effector target to the elbow-down joint-space target.
    x = target(1);
    y = target(2);
    c2 = (x^2 + y^2 - l1^2 - l2^2) / (2*l1*l2);
    c2 = max(min(c2, 1), -1);
    theta2 = acos(c2);
    k1 = l1 + l2*cos(theta2);
    k2 = l2*sin(theta2);
    theta1 = atan2(y, x) - atan2(k2, k1);
    xTargetJoint = [theta1; theta2; 0; 0];
end

function phi = encodeJointError(eJoint)
% Periodic joint-error encoding used by the network.
    phi = [sin(eJoint(1)); cos(eJoint(1)); sin(eJoint(2)); cos(eJoint(2)); eJoint(3); eJoint(4)];
end

function [xSim, finishTime] = simulateLqrArm(A, B, K, dt, t, xTargetJoint, target, l1, l2)
% Simulate the analytic LQR controller for the same horizon as the network.
    xSim = zeros(4, numel(t));
    consecutiveCounter = 0;
    finishTime = t(end);
    for k = 1:numel(t)-1
        eJoint = xSim(:, k) - xTargetJoint;
        u = -K * eJoint;
        xSim(:, k+1) = xSim(:, k) + dt*(A*xSim(:, k) + B*u);
        dist = endEffectorDistance(xSim(:, k+1), target, l1, l2);
        if dist < 0.1
            consecutiveCounter = consecutiveCounter + 1;
        else
            consecutiveCounter = 0;
        end
        if consecutiveCounter >= 10
            finishTime = (k+1)*dt;
            xSim(:, k+2:end) = repmat(xSim(:, k+1), 1, numel(t)-k-1);
            break
        end
    end
end

function [xSim, finishTime, finished] = simulateNetworkArm(net, A, B, dt, t, xTargetJoint, target, l1, l2, mu_X, sigma_X, mu_Y, sigma_Y, numInputs)
% Simulate the saved neural controller for a fresh reaching target.
    xSim = zeros(4, numel(t));
    consecutiveCounter = 0;
    finishTime = t(end);
    finished = false;
    for k = 1:numel(t)-1
        eJoint = xSim(:, k) - xTargetJoint;
        phi = encodeJointError(eJoint);
        netInput = (phi - mu_X) ./ (sigma_X * sqrt(numInputs));
        uNorm = numericOutput(predict(net, dlarray(netInput, 'CB')));
        u = uNorm .* sigma_Y + mu_Y;
        xSim(:, k+1) = xSim(:, k) + dt*(A*xSim(:, k) + B*u);
        dist = endEffectorDistance(xSim(:, k+1), target, l1, l2);
        if dist < 0.1
            consecutiveCounter = consecutiveCounter + 1;
        else
            consecutiveCounter = 0;
        end
        if consecutiveCounter >= 10
            finishTime = (k+1)*dt;
            xSim(:, k+2:end) = repmat(xSim(:, k+1), 1, numel(t)-k-1);
            finished = true;
            break
        end
    end
end

function dist = endEffectorDistance(xState, target, l1, l2)
% Distance from two-link end effector to target.
    th1 = xState(1);
    th2 = xState(2);
    endEffector = [l1*cos(th1) + l2*cos(th1 + th2); l1*sin(th1) + l2*sin(th1 + th2)];
    dist = hypot(endEffector(1) - target(1), endEffector(2) - target(2));
end

function M = localConcat(C, idx)
% Concatenate selected episode cells horizontally, skipping empties.
    M = [];
    for ii = 1:numel(idx)
        ci = C{idx(ii)};
        if ~isempty(ci)
            M = [M, ci]; %#ok<AGROW>
        end
    end
end

function assertNoSharedTargets(trainTargets, testTargets, message)
% Reject exact train/test target reuse after rounding for numerical stability.
    sharedTargets = intersect(round(trainTargets, 12), round(testTargets, 12), 'rows');
    assert(isempty(sharedTargets), message);
end

function assertNetworkIO(net, expectedInputs, expectedOutputs, label)
% Basic consistency check between the saved network and motor-control task.
    checkedLayers = false;
    try
        inputSize = net.Layers(1).InputSize;
        outputSize = net.Layers(end).OutputSize;
        checkedLayers = true;
    catch
    end
    if checkedLayers
        assert(inputSize(1) == expectedInputs, '%s input size mismatch.', label);
        assert(isequal(outputSize, expectedOutputs), '%s output size mismatch.', label);
    else
        y = numericOutput(predict(net, dlarray(zeros(expectedInputs, 1, 'single'), 'CB')));
        assert(numel(y) == expectedOutputs, '%s output size mismatch.', label);
    end
end`);

makeMlx("test_dynamical_system_saved_networks.mlx", `%% Test saved dynamical-system networks
% Loads each saved dynamical-system network with a dynamics-list entry and
% evaluates a perturbed closed-loop test trajectory. No training trajectory is
% used for metrics.

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
    net = loadSavedNetwork(modelPath);
    dynamicsName = erase(modelFiles(i), '_network.mat');
    rmse = dynamicalSystemTestRmse(repoRoot, net, dynamicsName);
    fprintf('%s perturbed test trajectory RMSE: %.4g\\n', dynamicsName, rmse);
end

${sharedHelpers}

function rmse = dynamicalSystemTestRmse(repoRoot, net, dynamicsName)
% Evaluate a saved network on a perturbed trajectory separate from training.
    dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
    names = string(dyn{:, 1});
    row = find(names == string(dynamicsName), 1);
    assert(~isempty(row), 'No dynamics-list row found for %s.', dynamicsName);
    d = dyn{row, 2};
    ic = dyn{row, 4:6}.';
    perturbationScale = 0.01;
    rng(9000 + row, 'twister');
    icTest = ic + perturbationScale * randn(size(ic));
    assert(norm(icTest - ic) > 0, '%s perturbation failed.', dynamicsName);
    tSpan = [0, 150];
    [~, xTestTrue] = ode45(@(t,y) int_dyn(y, string(dynamicsName), 0, 0, 0, 'Simulate'), ...
        tSpan, icTest(1:d));
    xTestTrue = normalize(xTestTrue) / sqrt(size(xTestTrue, 2));
    inputDelay = 1;
    xPred = xTestTrue(1:inputDelay, :);
    for k = 1:(size(xTestTrue, 1) - inputDelay)
        Xf = extractTimeDelayedInputsFeedback(xPred, inputDelay);
        yhat = predictDynamicalStep(net, Xf(:, end));
        xPred = [xPred; yhat(:).']; %#ok<AGROW>
    end
    assert(size(xPred, 1) == size(xTestTrue, 1), '%s prediction length mismatch.', dynamicsName);
    err = xPred - xTestTrue;
    rmse = sqrt(mean(err(:).^2, 'omitnan'));
end

function yhat = predictDynamicalStep(net, inCol)
% Run one closed-loop prediction step for a dynamical-system network.
    try
        yhat = predict(net, inCol(:).');
    catch
        yhat = predictCB(net, inCol(:));
    end
    yhat = numericOutput(yhat);
end`);

makeMlx("test_mnist_saved_networks.mlx", `%% Test saved MNIST-family network artifacts
% Loads each saved MNIST network artifact and evaluates only the dataset's
% official test split. No training or validation data are used for metrics.

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
    assert(isfile(networkPath), 'Missing saved network artifact: %s', networkPath);
    S = load(networkPath);
    if isfield(S, 'acc')
        fprintf('%s saved accuracy: %.2f%%\\n', networkFiles(i), accuracyPercent(S.acc));
    else
        fprintf('%s loaded; no acc variable stored.\\n', networkFiles(i));
    end
    net = loadSavedNetwork(networkPath);
    dataFile = erase(networkFiles(i), '_network.mat') + ".mat";
    dataPath = fullfile(repoRoot, 'data', dataFile);
    assert(isfile(dataPath), 'Missing MNIST test dataset: %s', dataPath);
    D = load(dataPath);
    assert(isfield(D, 'training') && isfield(D, 'test'), ...
        '%s must contain training and test structs.', dataFile);
    [XTest, YTest] = mnistTestData(D.training, D.test, dataFile);
    assertImageInputCompatible(net, XTest, dataFile);
    scores = numericOutput(predictImageBatch(net, XTest, 256));
    predIdx = scoresToClassIndex(scores, numel(categories(YTest)), numel(YTest));
    yHat = categorical(predIdx(:), 1:numel(categories(YTest)), categories(YTest));
    pct = 100 * mean(yHat(:) == YTest(:));
    fprintf('%s test-split accuracy: %.2f%%\\n', networkFiles(i), pct);
end

${sharedHelpers}

function [XTest, YTest] = mnistTestData(training, test, label)
% Load only the official test split and assert it is separate from training.
    assert(isfield(training, 'images') && isfield(training, 'labels'), ...
        '%s training split must contain images and labels.', label);
    assert(isfield(test, 'images') && isfield(test, 'labels'), ...
        '%s test split must contain images and labels.', label);
    nTrain = numel(training.labels);
    nTest = numel(test.labels);
    assert(nTrain > 0 && nTest > 0, '%s train/test split is empty.', label);
    assert(nTrain ~= nTest || ~isequal(size(training.images), size(test.images)) || ...
        ~isequal(training.images, test.images), ...
        '%s training and test image arrays are identical.', label);
    assertNoSharedImages(training.images, test.images, label);
    XTest = ensure4d(test.images);
    XTest = im2single(XTest);
    XTest = normalize(XTest, 4) / sqrt(28*28);
    XTest(isnan(XTest)) = 0;
    YTest = categorical(double(test.labels(:)) + 1);
end

function assertNoSharedImages(trainImages, testImages, label)
% Hard check for exact duplicate images across MNIST train/test splits.
    trainImages = ensure4d(trainImages);
    testImages = ensure4d(testImages);
    trainRows = reshape(trainImages, [], size(trainImages, 4)).';
    testRows = reshape(testImages, [], size(testImages, 4)).';
    sharedRows = intersect(trainRows, testRows, 'rows');
    assert(isempty(sharedRows), '%s exact duplicate images found across train/test.', label);
end

function assertImageInputCompatible(net, X, label)
% Fail before prediction if the test images do not match the saved network input.
    try
        inputSize = net.Layers(1).InputSize;
    catch
        return
    end
    if numel(inputSize) >= 3
        expected = inputSize(1:3);
        observed = [size(X, 1), size(X, 2), size(X, 3)];
        assert(isequal(observed, expected), ...
            '%s test image size %s does not match network input size %s.', ...
            label, mat2str(observed), mat2str(expected));
    end
end

function X = ensure4d(X)
% Normalize MNIST-family image arrays to H-by-W-by-C-by-N.
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

function Y = predictImageBatch(net, X, batchSize)
% Predict image batches while preserving the dataset test split.
    nObs = size(X, 4);
    Y = [];
    for first = 1:batchSize:nObs
        last = min(first + batchSize - 1, nObs);
        XBatch = X(:, :, :, first:last);
        try
            YBatch = predict(net, XBatch);
        catch
            YBatch = predict(net, dlarray(single(XBatch), 'SSCB'));
        end
        YBatch = numericOutput(YBatch);
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
% Convert class scores to class indices regardless of score orientation.
    if size(scores, 1) == nObs && size(scores, 2) == nClasses
        [~, idx] = max(scores, [], 2);
    elseif size(scores, 1) == nClasses && size(scores, 2) == nObs
        [~, idx] = max(scores, [], 1);
        idx = idx(:);
    else
        error('Unexpected score size %s for %d observations and %d classes.', ...
            mat2str(size(scores)), nObs, nClasses);
    end
end`);

makeMlx("test_pong_saved_network.mlx", `%% Test saved Pong network against a mirroring opponent
% This test uses the project Pong simulation test: the left paddle mirrors the
% ball position and the saved network controls the right paddle. The test loads
% the saved network from trained_networks and does not use training data.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

net = loadSavedNetwork(fullfile(repoRoot, 'trained_networks', 'Pong_network.mat'));
rng(1, 'twister');

% Internal truncated networks used for optional activation traces.
net2 = dlnetwork(net.Layers(1:5));
net3 = dlnetwork(net.Layers(1:3));
idx_neuron = [14068, 3042, 8509, 15830, 921]; %#ok<NASGU>

%% Simulation parameters
ballSpeed = 0.02;
maxBounceAngle = pi/4;
dt = 0.01;

paddleHeight = 0.2;
paddleWidth = 0.02;
ballRadius = 0.02;
oppPaddleX = 0.05;
netPaddleX = 1 - oppPaddleX - paddleWidth;

%% Scoring variables
mirrorScore = 0;
netContactCount = 0;
netMissCount = 0;

%% Simulation settings
maxNetworkHits = 100;
maxRounds = 500;
round = 0;

%% Set up animation
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

%% Main simulation loop
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
    set(oppPaddleRect, 'Position', [oppPaddleX, oppPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
    set(netPaddleRect, 'Position', [netPaddleX, netPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
    set(ballCircle, 'Position', [ballPos(1)-ballRadius, ballPos(2)-ballRadius, 2*ballRadius, 2*ballRadius]);
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
        pred = predict(net, state);
        pred = numericOutput(pred);

        [~, predClass] = max(pred, [], 2);
        if predClass == 1
            netAction = "Up";
        elseif predClass == 2
            netAction = "Down";
        else
            netAction = "Stay";
        end

        if netAction == "Up"
            netPaddleY = max(netPaddleY - 0.025, paddleHeight/2);
        elseif netAction == "Down"
            netPaddleY = min(netPaddleY + 0.025, 1 - paddleHeight/2);
        end

        ballPos = ballPos + ballVel;
        if ballPos(2) - ballRadius <= 0 || ballPos(2) + ballRadius >= 1
            ballVel(2) = -ballVel(2);
        end

        pred2 = numericOutput(predict(net2, state));
        pred1 = numericOutput(predict(net3, state));
        outy(cnt1) = netPaddleY; %#ok<SAGROW>
        truey(cnt1, :) = ballPos; %#ok<SAGROW>
        A2(cnt1, :) = pred2; %#ok<SAGROW>
        A1(cnt1, :) = pred1; %#ok<SAGROW>
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

        set(oppPaddleRect, 'Position', [oppPaddleX, oppPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
        set(netPaddleRect, 'Position', [netPaddleX, netPaddleY - paddleHeight/2, paddleWidth, paddleHeight]);
        set(ballCircle, 'Position', [ballPos(1)-ballRadius, ballPos(2)-ballRadius, 2*ballRadius, 2*ballRadius]);
        updateScore(scoreText, round, mirrorScore, netContactCount, netMissCount);
        drawnow;
        pause(dt);

        if netContactCount >= maxNetworkHits
            break
        end
    end
end

assert(netContactCount > 0 || netMissCount > 0, 'Pong simulation did not reach the network paddle.');
fprintf('After %d rounds:\n', round);
fprintf('Mirror (Left) Score: %d\n', mirrorScore);
fprintf('Network Hits: %d\n', netContactCount);
fprintf('Network Misses: %d\n', netMissCount);

${sharedHelpers}

function updateScore(scoreText, round, mirrorScore, netContactCount, netMissCount)
% Update the on-screen Pong scoreboard.
    scoreStr = sprintf('Round %d    Mirror (Left): %d    Network Hits: %d    Network Misses: %d', ...
        round, mirrorScore, netContactCount, netMissCount);
    set(scoreText, 'String', scoreStr);
end`);

console.log("Created MATLAB live-script tests in tests/.");
