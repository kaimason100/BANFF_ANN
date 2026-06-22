function run_dynamical_system_state_random_derivative_ode45(taskRow, weightSeed, trainingOptions)
%RUN_DYNAMICAL_SYSTEM_STATE_RANDOM_DERIVATIVE_ODE45 Train one DS derivative-field model.
% This variant keeps the state-random sampling, but trains the network in
% the older Dale-style form:
%   network(x) ~= x + dx/dt
% and validates with a fixed random derivative-field validation batch.
% Closed-loop ode45 rollouts are intentionally left for final testing.
% Outputs are saved to a distinct network set so previous results are untouched.
if nargin < 1 || isempty(taskRow)
    taskRow = 5;
end
if nargin < 2 || isempty(weightSeed)
    weightSeed = 0;
end
if nargin < 3 || isempty(trainingOptions)
    trainingOptions = struct();
end

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

totalEpochs = getOption(trainingOptions, 'TotalEpochs', 100000);
validationEveryEpochs = getOption(trainingOptions, 'ValidationEveryEpochs', 200);
validationSampleCount = getOption(trainingOptions, 'ValidationSamples', 100);
validationRandomSeed = getOption(trainingOptions, 'ValidationRandomSeed', 9000);
scalingSampleTime = getOption(trainingOptions, 'ScalingSampleTime', 1000);
randomSamplesPerEpoch = getOption(trainingOptions, 'RandomSamplesPerEpoch', 1000);
sampleRandomSeed = getOption(trainingOptions, 'SampleRandomSeed', 12345);
stateRandomRange = getOption(trainingOptions, 'StateRandomRange', [-1 1]);
saveProgressFigures = getOption(trainingOptions, 'SaveProgressFigures', true);
progressEveryEpochs = getOption(trainingOptions, 'ProgressEveryEpochs', 100);
learnRate = getOption(trainingOptions, 'LearnRate', 0.01);
gradientDecayFactor = getOption(trainingOptions, 'AdamGradientDecayFactor', 0.99);
squaredGradientDecayFactor = getOption(trainingOptions, 'AdamSquaredGradientDecayFactor', 0.999);
epsilon = getOption(trainingOptions, 'AdamEpsilon', 5e-6);
sourceNetworkSet = getOption(trainingOptions, 'SourceNetworkSet', '');
outputNetworkSet = getOption(trainingOptions, 'OutputNetworkSet', 'seeded_state_random_derivative_ode45_lr_0p01');
continuationLabel = getOption(trainingOptions, 'ContinuationLabel', '');
continueFromSource = getOption(trainingOptions, 'ContinueFromSource', false) || ~isempty(sourceNetworkSet);

trainDynamicalSystemDerivativeSeed(repoRoot, taskRow, weightSeed, totalEpochs, validationEveryEpochs, validationSampleCount, validationRandomSeed, saveProgressFigures, progressEveryEpochs, scalingSampleTime, randomSamplesPerEpoch, sampleRandomSeed, stateRandomRange, learnRate, gradientDecayFactor, squaredGradientDecayFactor, epsilon, continueFromSource, sourceNetworkSet, outputNetworkSet, continuationLabel);
end

function trainDynamicalSystemDerivativeSeed(repoRoot, ii, weightSeed, totalEpochs, validationEveryEpochs, validationSampleCount, validationRandomSeed, saveProgressFigures, progressEveryEpochs, scalingSampleTime, randomSamplesPerEpoch, sampleRandomSeed, stateRandomRange, learnRate, gradientDecayFactor, squaredGradientDecayFactor, epsilon, continueFromSource, sourceNetworkSet, outputNetworkSet, continuationLabel)
    weightCleanup = setWeightSeed(weightSeed); %#ok<NASGU>
    close all; clc; rng(1, 'twister');

    dyn = readtable(fullfile(repoRoot, 'examples', 'dynamical_systems', 'dynamics_list.xlsx'));
    dynamicsName = string(dyn{ii, 1});
    ic = dyn{ii, 4:6}.';
    d = dyn{ii, 2};
    rescaleFlag = dyn{ii, 8};
    assert(isfinite(scalingSampleTime) && scalingSampleTime > 0, '%s seed %d ScalingSampleTime must be positive.', dynamicsName, weightSeed);

    fprintf('\n[%s seed %d] preparing ode45 scaling sample over [0 %g]\n', ...
        dynamicsName, weightSeed, scalingSampleTime);
    [tTrain, xScaleRaw] = simulateReferenceScalingTrajectory(dynamicsName, ic(1:d), scalingSampleTime);
    trainStateRows = 1:size(xScaleRaw, 1);

    stateStats = fitReferenceStateScaling(xScaleRaw, rescaleFlag);
    assertReferenceStateScaling(stateStats, xScaleRaw, rescaleFlag, dynamicsName);
    xTrue = applyStateStandardization(xScaleRaw, stateStats);

    figure('Color', 'w');
    subplot(2,1,1); plot(tTrain, xTrue, 'LineWidth', 1); title(sprintf('%s normalized true time series', dynamicsName)); xlabel('Time'); ylabel('State');
    if size(xTrue, 2) > 1
        subplot(2,1,2); plot(xTrue(:,1), xTrue(:,2), 'LineWidth', 1); title('Training/validation phase plot'); xlabel('x_1'); ylabel('x_2');
    end

    randomSamplesPerEpoch = double(randomSamplesPerEpoch);
    assert(isfinite(randomSamplesPerEpoch) && randomSamplesPerEpoch == round(randomSamplesPerEpoch) && randomSamplesPerEpoch >= 1, ...
        '%s seed %d RandomSamplesPerEpoch must be a positive integer.', dynamicsName, weightSeed);
    assert(numel(stateRandomRange) == 2 && all(isfinite(stateRandomRange)) && stateRandomRange(2) > stateRandomRange(1), ...
        '%s seed %d StateRandomRange must be [min max].', dynamicsName, weightSeed);
    initialSampleRng = RandStream('mt19937ar', 'Seed', double(sampleRandomSeed));
    [XTrain, YTrain] = randomStateSpaceDerivativeTrainingSamples(dynamicsName, stateStats, d, randomSamplesPerEpoch, initialSampleRng, stateRandomRange);
    validationRng = RandStream('mt19937ar', 'Seed', double(validationRandomSeed));
    [XValidation, YValidation] = randomStateSpaceDerivativeTrainingSamples(dynamicsName, stateStats, d, double(validationSampleCount), validationRng, stateRandomRange);
    fprintf('[%s seed %d] fixed state-random derivative-field training | scaled range [%g %g] | training samples: %d | training RNG seed: %d\n', ...
        dynamicsName, weightSeed, stateRandomRange(1), stateRandomRange(2), randomSamplesPerEpoch, double(sampleRandomSeed));
    fprintf('[%s seed %d] derivative-field validation | scaled range [%g %g] | fixed validation samples: %d | validation RNG seed: %d\n', ...
        dynamicsName, weightSeed, stateRandomRange(1), stateRandomRange(2), double(validationSampleCount), double(validationRandomSeed));
    fprintf('[%s seed %d] reference scaling | rescale flag=%d | mean=[%s] | scale diag=[%s]\n', ...
        dynamicsName, weightSeed, double(rescaleFlag), num2str(stateStats.mu, ' %.6g'), num2str(stateStats.sigma, ' %.6g'));

    layers = [
        featureInputLayer(size(XTrain, 2), Normalization='none')
        fullyConnectedLayer(16000, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
        tanhLayer
        fullyConnectedLayer(16000, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
        tanhLayer
        fullyConnectedLayer(size(xTrue, 2), 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')];

    validationEveryEpochs = max(1, min(validationEveryEpochs, totalEpochs));
    nChecks = ceil(totalEpochs / validationEveryEpochs) + double(continueFromSource);
    validationHistory = table('Size', [nChecks, 4], ...
        'VariableTypes', {'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Epoch', 'ValidationLoss', 'TrainingLoss', 'NumValidationSamples'});
    safeDynamicsName = matlab.lang.makeValidName(char(dynamicsName));
    safeOutputSet = matlab.lang.makeValidName(char(outputNetworkSet));
    progressFigureDir = fullfile(repoRoot, 'outputs', 'state_random_derivative_ode45_validation_progress', safeOutputSet, safeDynamicsName, sprintf('seed_%03d', double(weightSeed)));
    progressCheckpointDir = fullfile(repoRoot, 'outputs', 'state_random_derivative_ode45_validation_checkpoints', safeOutputSet, safeDynamicsName, sprintf('seed_%03d', double(weightSeed)));
    if exist(progressFigureDir, 'dir') ~= 7, mkdir(progressFigureDir); end
    if exist(progressCheckpointDir, 'dir') ~= 7, mkdir(progressCheckpointDir); end

    [useGPU, executionDevice] = selectExecutionDevice();
    fprintf('[%s seed %d] training execution device: %s\n', dynamicsName, weightSeed, executionDevice);

    sourcePath = '';
    sourceMetadata = struct();
    sourceValidationInfo = struct();
    if continueFromSource
        assert(~isempty(sourceNetworkSet), '%s seed %d continuation requested without SourceNetworkSet.', dynamicsName, weightSeed);
        sourcePath = stateRandomDerivativeSeededNetworkPath(repoRoot, char(dynamicsName), weightSeed, sourceNetworkSet, false);
        assert(isfile(sourcePath), 'Missing source derivative-field network for continuation: %s', sourcePath);
        sourceSave = load(sourcePath);
        assert(isfield(sourceSave, 'net') && isa(sourceSave.net, 'dlnetwork'), '%s seed %d source file does not contain dlnetwork net.', dynamicsName, weightSeed);
        assert(isfield(sourceSave, 'stateStats'), '%s seed %d source file lacks stateStats.', dynamicsName, weightSeed);
        assertStatsClose(sourceSave.stateStats.mu, stateStats.mu, dynamicsName, 'source state mean');
        assertStatsClose(sourceSave.stateStats.sigma, stateStats.sigma, dynamicsName, 'source state max-absolute scale');
        assertStatsClose(sourceSave.stateStats.scale, stateStats.scale, dynamicsName, 'source state scale');
        currentNet = sourceSave.net;
        if isfield(sourceSave, 'metadata'), sourceMetadata = sourceSave.metadata; end
        if isfield(sourceSave, 'validationInfo'), sourceValidationInfo = sourceSave.validationInfo; end
        fprintf('[%s seed %d] continuing from source network set %s: %s\n', dynamicsName, weightSeed, char(sourceNetworkSet), sourcePath);
    else
        currentNet = dlnetwork(layerGraph(layers));
    end
    if useGPU
        currentNet = moveDlnetwork(currentNet, @gpuArray);
    end
    biasMask = strcmp(currentNet.Learnables.Parameter, 'Bias');
    biasValues = currentNet.Learnables.Value(biasMask);
    bestNet = [];
    bestValidationLoss = Inf;
    bestEpoch = NaN;
    trainingInfo = table('Size', [totalEpochs, 3], ...
        'VariableTypes', {'double', 'double', 'double'}, ...
        'VariableNames', {'Epoch', 'TrainingLoss', 'NumSamples'});
    trailingAvg = [];
    trailingAvgSq = [];
    iteration = 0;
    checkIdx = 0;
    progressTimer = tic;

    lossFig = figure('Color', 'w');
    lossAx = axes(lossFig);
    lossLine = animatedline(lossAx, 'Color', [0 0.4470 0.7410], 'LineWidth', 1);
    xlabel(lossAx, 'Epoch'); ylabel(lossAx, 'Training loss'); grid(lossAx, 'on');
    set(lossAx, 'YScale', 'log');
    title(lossAx, sprintf('%s seed %d continuous Adam training loss', dynamicsName, double(weightSeed)));

    validationLossPlot = initializeLiveValidationLossPlot(dynamicsName, weightSeed);
    if useGPU
        dlXTrain = dlarray(gpuArray(single(XTrain.')), 'CB');
        dlYTrain = dlarray(gpuArray(single(YTrain.')), 'CB');
        dlXValidation = dlarray(gpuArray(single(XValidation.')), 'CB');
        dlYValidation = dlarray(gpuArray(single(YValidation.')), 'CB');
    else
        dlXTrain = dlarray(single(XTrain.'), 'CB');
        dlYTrain = dlarray(single(YTrain.'), 'CB');
        dlXValidation = dlarray(single(XValidation.'), 'CB');
        dlYValidation = dlarray(single(YValidation.'), 'CB');
    end

    if continueFromSource
        checkIdx = checkIdx + 1;
        sourceValidationLoss = derivativeFieldLoss(currentNet, dlXValidation, dlYValidation);
        sourceValidationLossScalar = double(gather(extractdata(sourceValidationLoss)));
        updateLiveValidationLossPlot(validationLossPlot, 0, sourceValidationLossScalar);
        validationHistory.Epoch(checkIdx) = 0;
        validationHistory.ValidationLoss(checkIdx) = sourceValidationLossScalar;
        validationHistory.TrainingLoss(checkIdx) = NaN;
        validationHistory.NumValidationSamples(checkIdx) = validationSampleCount;
        if isfinite(sourceValidationLossScalar)
            bestValidationLoss = sourceValidationLossScalar;
            bestEpoch = 0;
            bestNet = currentNet;
        end
        fprintf('[%s seed %d] source epoch 0 validation loss %.6g\n', dynamicsName, weightSeed, sourceValidationLossScalar);
    end

    for epoch = 1:totalEpochs
        iteration = iteration + 1;
        [gradients, lossValue] = dlfeval(@modelGradients, currentNet, dlXTrain, dlYTrain, biasMask);
        [biasValues, trailingAvg, trailingAvgSq] = adamUpdateBiasValues(biasValues, gradients, trailingAvg, trailingAvgSq, ...
            iteration, learnRate, gradientDecayFactor, squaredGradientDecayFactor, epsilon);
        currentNet = setBiasLearnables(currentNet, biasMask, biasValues);

        lossScalar = double(gather(extractdata(lossValue)));
        trainingInfo.Epoch(epoch) = epoch;
        trainingInfo.TrainingLoss(epoch) = lossScalar;
        trainingInfo.NumSamples(epoch) = randomSamplesPerEpoch;
        addpoints(lossLine, epoch, positiveForLogPlot(lossScalar));
        if mod(epoch, 10) == 0 || epoch == 1 || epoch == totalEpochs
            drawnow limitrate;
        end
        if mod(epoch, progressEveryEpochs) == 0 || epoch == 1 || epoch == totalEpochs
            printDerivativeTrainingProgress(dynamicsName, weightSeed, epoch, totalEpochs, lossScalar, progressTimer);
        end

        if mod(epoch, validationEveryEpochs) == 0 || epoch == totalEpochs
            checkIdx = checkIdx + 1;
            validationLoss = derivativeFieldLoss(currentNet, dlXValidation, dlYValidation);
            validationLossScalar = double(gather(extractdata(validationLoss)));
            updateLiveValidationLossPlot(validationLossPlot, epoch, validationLossScalar);
            validationHistory.Epoch(checkIdx) = epoch;
            validationHistory.ValidationLoss(checkIdx) = validationLossScalar;
            validationHistory.TrainingLoss(checkIdx) = lossScalar;
            validationHistory.NumValidationSamples(checkIdx) = validationSampleCount;
            if saveProgressFigures
                saveValidationProgressOutputs(lossFig, validationLossPlot.fig, progressFigureDir, progressCheckpointDir, dynamicsName, weightSeed, epoch, trainingInfo, validationHistory, checkIdx);
            end

            if isfinite(validationLossScalar) && validationLossScalar < bestValidationLoss
                bestValidationLoss = validationLossScalar;
                bestEpoch = epoch;
                bestNet = currentNet;
            end
            fprintf('[%s seed %d] epoch %d/%d | validation loss %.6g | best validation loss %.6g at epoch %d\n', ...
                dynamicsName, weightSeed, epoch, totalEpochs, validationLossScalar, bestValidationLoss, bestEpoch);
        end
    end

    assert(~isempty(bestNet), '%s seed %d never produced a finite derivative-field validation loss.', dynamicsName, weightSeed);
    net = moveDlnetwork(bestNet, @gather);
    bestHistoryRow = find(validationHistory.Epoch == bestEpoch, 1, 'first');
    assert(~isempty(bestHistoryRow), '%s seed %d best epoch was not found in validation history.', dynamicsName, weightSeed);
    validationInfo.bestEpoch = bestEpoch;
    validationInfo.bestValidationLoss = bestValidationLoss;
    validationInfo.bestSelectionMetric = bestValidationLoss;
    validationInfo.history = validationHistory;
    validationInfo.selectionMetric = 'DerivativeFieldValidationLoss';
    validationInfo.selectionAggregation = 'mean';
    validationInfo.selectionRule = 'minimum fixed random derivative-field validation mean squared error';
    validationInfo.validationEveryEpochs = validationEveryEpochs;
    validationInfo.totalEpochs = totalEpochs;
    validationInfo.validationRandomSeed = validationRandomSeed;
    validationInfo.validationSampleCount = validationSampleCount;
    validationInfo.validationSampleRandomSeed = double(validationRandomSeed);
    validationInfo.progressFigureDir = progressFigureDir;
    validationInfo.progressCheckpointDir = progressCheckpointDir;
    validationInfo.scalingSampleTime = scalingSampleTime;
    validationInfo.scalingRescaleFlag = rescaleFlag;
    validationInfo.randomSamplesPerEpoch = randomSamplesPerEpoch;
    validationInfo.sampleRandomSeed = sampleRandomSeed;
    validationInfo.fixedTrainingSamples = true;
    validationInfo.sameTrainingValidationDataAcrossSeeds = true;
    validationInfo.trainingSamplingMode = 'state-random-derivative-field';
    validationInfo.derivativeTargetForm = 'network(x_norm) = x_norm + dx_norm_dt';
    validationInfo.stateRandomRange = stateRandomRange;
    validationInfo.learningRate = learnRate;
    validationInfo.adamGradientDecayFactor = gradientDecayFactor;
    validationInfo.adamSquaredGradientDecayFactor = squaredGradientDecayFactor;
    validationInfo.adamEpsilon = epsilon;
    validationInfo.integrationMethod = 'ode45-learned-derivative-field';
    validationInfo.trueTrajectoryIntegrationMethod = 'ode45-reference-scaling-sample';
    validationInfo.continuation = continueFromSource;
    validationInfo.sourceNetworkSet = char(sourceNetworkSet);
    validationInfo.outputNetworkSet = char(outputNetworkSet);
    validationInfo.sourceNetworkPath = char(sourcePath);
    validationInfo.continuationLabel = char(continuationLabel);
    validationInfo.sourceValidationInfo = sourceValidationInfo;

    split.trainTimeRows = trainStateRows;
    split.validation = 'fixed random reference-scaled state-space derivative-field samples';
    metadata.task = char(dynamicsName);
    metadata.normalization = 'reference mean/max-absolute trajectory scaling';
    metadata.seed = weightSeed;
    metadata.seedList = 0:9;
    metadata.selection = 'derivative-field-validation-loss';
    metadata.selectionMetric = 'DerivativeFieldValidationLoss';
    metadata.trainingImplementation = 'continuous custom Adam loop with fixed uniformly random reference-scaled state-space samples, Dale-style derivative-field targets, and fixed random validation loss selection';
    metadata.rngSeed = 1;
    metadata.validationRandomSeed = validationRandomSeed;
    metadata.validationSampleCount = validationSampleCount;
    metadata.validationSampleRandomSeed = double(validationRandomSeed);
    metadata.scalingSampleTime = scalingSampleTime;
    metadata.scalingRescaleFlag = rescaleFlag;
    metadata.randomSamplesPerEpoch = randomSamplesPerEpoch;
    metadata.sampleRandomSeed = sampleRandomSeed;
    metadata.fixedTrainingSamples = true;
    metadata.sameTrainingValidationDataAcrossSeeds = true;
    metadata.trainingSamplingMode = 'state-random-derivative-field';
    metadata.derivativeTargetForm = 'network(x_norm) = x_norm + dx_norm_dt';
    metadata.stateRandomRange = stateRandomRange;
    metadata.learningRate = learnRate;
    metadata.adamGradientDecayFactor = gradientDecayFactor;
    metadata.adamSquaredGradientDecayFactor = squaredGradientDecayFactor;
    metadata.adamEpsilon = epsilon;
    metadata.integrationMethod = 'ode45-learned-derivative-field';
    metadata.trueTrajectoryIntegrationMethod = 'ode45-reference-scaling-sample';
    metadata.networkSet = char(outputNetworkSet);
    metadata.progressEveryEpochs = progressEveryEpochs;
    metadata.saveProgressFigures = saveProgressFigures;
    metadata.executionDevice = executionDevice;
    metadata.arcCompatible = true;
    metadata.continuation = continueFromSource;
    metadata.sourceNetworkSet = char(sourceNetworkSet);
    metadata.sourceNetworkPath = char(sourcePath);
    metadata.continuationLabel = char(continuationLabel);
    metadata.sourceMetadata = sourceMetadata;

    modelPath = stateRandomDerivativeSeededNetworkPath(repoRoot, char(dynamicsName), weightSeed, outputNetworkSet, true);
    save(modelPath, 'net', 'trainingInfo', 'split', 'stateStats', 'metadata', 'validationInfo', '-v7.3');
    fprintf('[%s seed %d] saved best derivative-field model from epoch %d to %s\n', dynamicsName, weightSeed, bestEpoch, modelPath);
end


function value = getOption(options, fieldName, defaultValue)
    value = defaultValue;
    if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    end
end

function printDerivativeTrainingProgress(dynamicsName, weightSeed, epochValue, totalEpochs, lossScalar, progressTimer)
    elapsedSeconds = toc(progressTimer);
    pct = 100 * double(epochValue) / double(totalEpochs);
    if epochValue > 0
        etaSeconds = elapsedSeconds * (double(totalEpochs) / double(epochValue) - 1);
    else
        etaSeconds = NaN;
    end
    fprintf('[derivative-field %s seed %d progress] epoch %d/%d | %.1f%% complete | train loss %.6g | elapsed %s | ETA %s\n', ...
        dynamicsName, double(weightSeed), double(epochValue), double(totalEpochs), pct, lossScalar, formatDurationSeconds(elapsedSeconds), formatDurationSeconds(etaSeconds));
end

function txt = formatDurationSeconds(secondsValue)
    if ~isfinite(secondsValue) || secondsValue < 0
        txt = 'n/a';
        return
    end
    secondsValue = round(secondsValue);
    hoursValue = floor(secondsValue / 3600);
    minutesValue = floor(mod(secondsValue, 3600) / 60);
    secsValue = mod(secondsValue, 60);
    txt = sprintf('%02d:%02d:%02d', hoursValue, minutesValue, secsValue);
end

function saveValidationProgressOutputs(lossFig, validationFig, progressFigureDir, progressCheckpointDir, dynamicsName, weightSeed, epochValue, trainingInfo, validationHistory, checkIdx)
    safeName = matlab.lang.makeValidName(char(dynamicsName));
    lossLatest = fullfile(progressFigureDir, sprintf('%s_seed_%03d_loss_latest.png', safeName, double(weightSeed)));
    validationLatest = fullfile(progressFigureDir, sprintf('%s_seed_%03d_validation_latest.png', safeName, double(weightSeed)));
    lossEpoch = fullfile(progressFigureDir, sprintf('%s_seed_%03d_epoch_%06d_loss.png', safeName, double(weightSeed), double(epochValue)));
    validationEpoch = fullfile(progressFigureDir, sprintf('%s_seed_%03d_epoch_%06d_validation.png', safeName, double(weightSeed), double(epochValue)));
    saveas(lossFig, lossLatest);
    saveas(validationFig, validationLatest);
    saveas(lossFig, lossEpoch);
    saveas(validationFig, validationEpoch);

    trainingInfoPartial = trainingInfo(1:epochValue, :); %#ok<NASGU>
    validationHistoryPartial = validationHistory(1:checkIdx, :); %#ok<NASGU>
    latestCheckpoint = fullfile(progressCheckpointDir, sprintf('%s_seed_%03d_latest_progress.mat', safeName, double(weightSeed)));
    epochCheckpoint = fullfile(progressCheckpointDir, sprintf('%s_seed_%03d_epoch_%06d_progress.mat', safeName, double(weightSeed), double(epochValue)));
    save(latestCheckpoint, 'trainingInfoPartial', 'validationHistoryPartial', '-v7.3');
    save(epochCheckpoint, 'trainingInfoPartial', 'validationHistoryPartial', '-v7.3');
end

function [useGPU, executionDevice] = selectExecutionDevice()
    useGPU = false;
    executionDevice = 'CPU';
    try
        useGPU = canUseGPU();
    catch
        try
            gpuDevice();
            useGPU = true;
        catch
            useGPU = false;
        end
    end
    if useGPU
        try
            g = gpuDevice();
            executionDevice = sprintf('GPU: %s', g.Name);
        catch
            executionDevice = 'GPU';
        end
    end
end

function net = moveDlnetwork(net, moveFcn)
    learnables = net.Learnables;
    for k = 1:numel(learnables.Value)
        v = learnables.Value{k};
        if isa(v, 'dlarray')
            learnables.Value{k} = dlarray(moveFcn(extractdata(v)));
        else
            learnables.Value{k} = moveFcn(v);
        end
    end
    net.Learnables = learnables;
end

function validationLossPlot = initializeLiveValidationLossPlot(dynamicsName, weightSeed)
    fig = figure('Color', 'w');
    ax = axes(fig);
    hold(ax, 'on'); grid(ax, 'on');
    validationLine = animatedline(ax, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
    xlabel(ax, 'Epoch'); ylabel(ax, 'Validation loss');
    set(ax, 'YScale', 'log');
    title(ax, sprintf('%s seed %d derivative-field validation loss', dynamicsName, double(weightSeed)));
    validationLossPlot.fig = fig;
    validationLossPlot.ax = ax;
    validationLossPlot.line = validationLine;
end

function updateLiveValidationLossPlot(validationLossPlot, epochValue, validationLoss)
    addpoints(validationLossPlot.line, epochValue, positiveForLogPlot(validationLoss));
    drawnow limitrate;
end

function y = positiveForLogPlot(x)
    y = double(x);
    y(~isfinite(y) | y <= 0) = realmin('double');
end

function [XBatch, YBatch] = randomStateSpaceDerivativeTrainingSamples(dynamicsName, stateStats, stateDim, nRequestedSamples, rngStream, stateRandomRange)
    assert(nRequestedSamples >= 1, '%s state-random training sample count must be positive.', dynamicsName);
    lo = double(stateRandomRange(1));
    hi = double(stateRandomRange(2));
    XBatch = lo + (hi - lo) * rand(rngStream, nRequestedSamples, stateDim);
    rawX = inverseStateStandardization(XBatch, stateStats);
    rawDX = dynamicsDerivative(dynamicsName, rawX.').';
    dXNormDt = stateStats.scale * (rawDX ./ stateStats.sigma);
    YBatch = XBatch + dXNormDt;
    assert(all(isfinite(YBatch(:))), '%s state-random derivative-field targets contain non-finite values.', dynamicsName);
end

function [gradients, lossValue] = modelGradients(net, dlX, dlY, biasMask)
    lossValue = derivativeFieldLoss(net, dlX, dlY);
    biasValues = net.Learnables.Value(biasMask);
    gradients = dlgradient(lossValue, biasValues);
end

function lossValue = derivativeFieldLoss(net, dlX, dlY)
    dlYPred = forward(net, dlX);
    lossValue = mean((dlYPred - dlY).^2, 'all');
end

function [biasValues, trailingAvg, trailingAvgSq] = adamUpdateBiasValues(biasValues, gradients, trailingAvg, trailingAvgSq, iteration, learnRate, beta1, beta2, epsilon)
    if isempty(trailingAvg)
        trailingAvg = cell(size(biasValues));
        trailingAvgSq = cell(size(biasValues));
        for k = 1:numel(biasValues)
            trailingAvg{k} = zeros(size(biasValues{k}), 'like', biasValues{k});
            trailingAvgSq{k} = zeros(size(biasValues{k}), 'like', biasValues{k});
        end
    end

    for k = 1:numel(biasValues)
        g = gradients{k};
        trailingAvg{k} = beta1 * trailingAvg{k} + (1 - beta1) * g;
        trailingAvgSq{k} = beta2 * trailingAvgSq{k} + (1 - beta2) * (g.^2);
        avgCorrected = trailingAvg{k} ./ (1 - beta1^iteration);
        avgSqCorrected = trailingAvgSq{k} ./ (1 - beta2^iteration);
        biasValues{k} = biasValues{k} - learnRate * avgCorrected ./ (sqrt(avgSqCorrected) + epsilon);
    end
end

function net = setBiasLearnables(net, biasMask, biasValues)
    learnables = net.Learnables;
    learnables.Value(biasMask) = biasValues;
    net.Learnables = learnables;
end

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
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'src')));
end

function [t, x] = simulateReferenceScalingTrajectory(dynamicsName, initialCondition, scalingSampleTime)
    [t, x] = ode45(@(tValue, xValue) dynamicsDerivative(dynamicsName, xValue), [0 double(scalingSampleTime)], double(initialCondition(:)));
    assert(all(isfinite(x(:))), '%s ode45 scaling sample contains non-finite values.', dynamicsName);
end

function stats = fitReferenceStateScaling(xSample, rescaleFlag)
    if double(rescaleFlag) == 1
        stats.sigma = max(abs(xSample), [], 1);
        stats.mu = mean(xSample, 1);
    else
        stats.sigma = ones(1, size(xSample, 2));
        stats.mu = zeros(1, size(xSample, 2));
    end
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1;
    stats.method = 'reference mean/max-absolute trajectory scaling';
end

function x = applyStateStandardization(x, stats)
    x = stats.scale * ((x - stats.mu) ./ stats.sigma);
    x(~isfinite(x)) = 0;
end

function x = inverseStateStandardization(xNorm, stats)
    x = (xNorm ./ stats.scale) .* stats.sigma + stats.mu;
end

function assertReferenceStateScaling(stats, xSample, rescaleFlag, label)
    expected = fitReferenceStateScaling(xSample, rescaleFlag);
    assertStatsClose(stats.mu, expected.mu, label, 'state mean');
    assertStatsClose(stats.sigma, expected.sigma, label, 'state max-absolute scale');
    assertStatsClose(stats.scale, expected.scale, label, 'state scale');
end

function assertStatsClose(actual, expected, label, statName)
    actual = double(actual);
    expected = double(expected);
    assert(isequal(size(actual), size(expected)), '%s %s size mismatch.', label, statName);
    tol = 1e-10 * max(1, max(abs(expected(:))));
    assert(all(abs(actual(:) - expected(:)) <= tol), '%s %s was not fitted from the training split.', label, statName);
end

function dx = dynamicsDerivative(dynamicsName, x)
    inputSize = size(x);
    if isvector(x)
        x = x(:);
        inputSize = size(x);
    end
    dx = int_dyn(x, string(dynamicsName), 0, 0, 0, 'Simulate');
    assert(isequal(size(dx), inputSize), '%s dynamics derivative returned size [%s] for input size [%s].', ...
        dynamicsName, num2str(size(dx)), num2str(inputSize));
end

function cleanupObj = setWeightSeed(seedValue)
    setappdata(0, 'NAR_WEIGHT_SEED_BASE', double(seedValue));
    fprintf('[seeded derivative-field training] NAR_WEIGHT_SEED_BASE=%d\n', double(seedValue));
    cleanupObj = onCleanup(@() clearWeightSeed());
end

function clearWeightSeed()
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        rmappdata(0, 'NAR_WEIGHT_SEED_BASE');
    end
end

function modelPath = stateRandomDerivativeSeededNetworkPath(repoRoot, taskName, seedValue, networkSet, createFolder)
    if nargin < 4 || isempty(networkSet)
        networkSet = 'seeded_state_random_derivative_ode45_lr_0p01';
    end
    if nargin < 5 || isempty(createFolder)
        createFolder = true;
    end
    folderName = matlab.lang.makeValidName(char(taskName));
    outDir = fullfile(repoRoot, 'trained_networks', char(networkSet), folderName);
    if createFolder && exist(outDir, 'dir') ~= 7, mkdir(outDir); end
    modelPath = fullfile(outDir, sprintf('%s_seed_%03d_network.mat', folderName, double(seedValue)));
end
