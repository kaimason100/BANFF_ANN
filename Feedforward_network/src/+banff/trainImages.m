function trainImages(taskName, SEEDS, cfg)
% Shared numerical implementation for local and Slurm entry points.
repoRoot = banff.root();
OUT_ROOT = fullfile(repoRoot,'trained_networks',cfg.NetworkSet);
DATA_DIR = fullfile(repoRoot,'data');
WIDTH = cfg.Width; NUM_EPOCHS = cfg.MaxEpochs;
LEARN_RATE = cfg.LearnRate; BATCH_SIZE = cfg.MiniBatchSize; SPLIT_SEED = cfg.SplitSeed;
targetList = string(taskName) + ".mat";
for ds = targetList
    ds = char(ds);
    dataPath = fullfile(DATA_DIR, ds);
    assert(isfile(dataPath), 'Missing MNIST-family dataset: %s', dataPath);
    D = load(dataPath);
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = banff.shared.mnistData(D.training, D.test, ds);
    outName = erase(ds, '.mat');
    for seedValue = SEEDS
        fprintf('\n[%s] Training seed %d\n', outName, seedValue);
        banff.shared.setWeightSeed(seedValue);
        try
        rng(SPLIT_SEED, 'twister');
        [trainIdx, valIdx] = banff.imageSplit(XTrainFull, YTrainFull, SPLIT_SEED);
        imageStats = banff.shared.fitImageStandardization(XTrainFull(:, :, :, trainIdx));
        banff.shared.assertImageStatsFromTraining(imageStats, XTrainFull, trainIdx, ds);
        banff.assertImagePartitions(XTrainFull,XTest,imageStats,trainIdx,valIdx);
        XTrain = banff.shared.applyImageStandardization(XTrainFull(:, :, :, trainIdx), imageStats); YTrain = YTrainFull(trainIdx);
        XVal = banff.shared.applyImageStandardization(XTrainFull(:, :, :, valIdx), imageStats); YVal = YTrainFull(valIdx);
        layers = [
            imageInputLayer([28 28 1], Normalization='none')
            fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights1_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
            tanhLayer
            fullyConnectedLayer(WIDTH, 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights2_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
            tanhLayer
            fullyConnectedLayer(numOut, 'Name', 'output', 'WeightLearnRateFactor', 0, 'WeightsInitializer', @customWeights3_git, 'BiasL2Factor', 0, 'BiasInitializer', 'zeros')
            softmaxLayer
    classificationLayer];
                options = trainingOptions('adam', ...
            'MaxEpochs', NUM_EPOCHS, 'MiniBatchSize', BATCH_SIZE, ...
            'InitialLearnRate', LEARN_RATE, 'Shuffle', 'every-epoch', ...
            'ValidationData', {XVal, YVal}, 'ValidationFrequency', 10, 'OutputNetwork', 'best-validation-loss', ...
            'Plots', banff.plotOption(cfg.Plots), 'ExecutionEnvironment', cfg.ExecutionEnvironment, 'Verbose', 0, ...
            'OutputFcn', @(info) banff.progress(info,cfg.MaxEpochs,char(taskName)));
        [net, info] = trainNetwork(XTrain, YTrain, layers, options);
        XTestEval = banff.shared.applyImageStandardization(XTest, imageStats);
        scores = banff.shared.numericOutput(banff.imageTraining.predictImageBatches(net, XTestEval, 256));
        predIdx = banff.shared.scoresToClassIndex(scores, numOut, numel(YTest));
        yHat = categorical(predIdx(:), 1:numOut);
        acc = 100 * mean(yHat(:) == YTest(:));
        split.trainIdx = trainIdx; split.valIdx = valIdx; split.testSplit = 'official_test_struct';
        metadata.task = outName; metadata.dataset = ds; metadata.seed = seedValue; metadata.seedList = 0:9; metadata.width = WIDTH;
        metadata.normalization = 'train-fitted image z-score';
        metadata.datasetSHA256 = banff.fileSHA256(dataPath);
        metadata.trainingOptions = cfg;
        outDir = fullfile(OUT_ROOT, outName); if ~exist(outDir, 'dir'), mkdir(outDir); end
        save(fullfile(outDir, sprintf('%s_seed_%03d_network.mat', outName, seedValue)), 'net', 'info', 'acc', 'split', 'imageStats', 'metadata', '-v7.3');
        banff.shared.clearWeightSeed();
        fprintf('[%s] seed %d test accuracy: %.2f%%\n', outName, seedValue, acc);
        catch ME
            banff.shared.clearWeightSeed();
            rethrow(ME);
        end
    end
end



end
