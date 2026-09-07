function testImages(taskName, SEEDS)
%% Aggregate seeded MNIST-family classification test results
repoRoot = banff.root();


DATA_DIR = fullfile(repoRoot, 'data');
targetList = string(taskName) + ".mat";

for ds = targetList
    ds = char(ds);
    outName = erase(ds, '.mat');
    D = load(fullfile(DATA_DIR, ds));
    [XTrainFull, YTrainFull, XTest, YTest, numOut] = banff.shared.mnistData(D.training, D.test, ds);
    [expectedTrainIdx, expectedValIdx] = banff.imageSplit(XTrainFull, YTrainFull, 42);
    acc = nan(numel(SEEDS), 1);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded_grouped', outName, sprintf('%s_seed_%03d_network.mat', outName, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assert(isfield(S.metadata,'datasetSHA256') && strcmp(S.metadata.datasetSHA256, ...
            banff.fileSHA256(fullfile(DATA_DIR,ds))), 'Dataset provenance mismatch; retrain with this branch.');
        banff.shared.assertSeedMetadata(S, outName, seedValue, SEEDS);
        assert(isfield(S, 'imageStats'), 'MNIST seeded model lacks train-fitted image normalization metadata.');
        assert(isfield(S, 'split') && strcmp(S.split.testSplit, 'official_test_struct'), 'MNIST seeded model lacks test-split metadata.');
        assert(isfield(S.split, 'trainIdx') && isfield(S.split, 'valIdx'), 'MNIST seeded model lacks train/validation split metadata.');
        assert(isempty(intersect(S.split.trainIdx(:), S.split.valIdx(:))), '%s seed %d train/validation leakage.', outName, seedValue);
        assert(all(S.split.trainIdx(:) >= 1 & S.split.trainIdx(:) <= numel(YTrainFull)) && all(S.split.valIdx(:) >= 1 & S.split.valIdx(:) <= numel(YTrainFull)), ...
            '%s seed %d train/validation split indices are out of range.', outName, seedValue);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)), ...
            '%s seed %d train/validation split metadata does not match the deterministic training split.', outName, seedValue);
        banff.shared.assertImageStatsFromTraining(S.imageStats, XTrainFull, S.split.trainIdx, outName);
        banff.assertImagePartitions(XTrainFull,XTest,S.imageStats,S.split.trainIdx,S.split.valIdx);
        net = banff.shared.loadSavedNetwork(modelPath);
        XTestEval = banff.shared.applyImageStandardization(XTest, S.imageStats);
        scores = banff.shared.numericOutput(banff.imageTesting.predictImageBatches(net, XTestEval, 256));
        predIdx = banff.shared.scoresToClassIndex(scores, numOut, numel(YTest));
        yHat = categorical(predIdx(:), 1:numOut);
        acc(i) = mean(yHat(:) == YTest(:));
        fprintf('%s seed %d accuracy: %.2f%%\n', outName, seedValue, 100*acc(i));
    end
    fprintf('%s seeded accuracy mean ± SD: %.2f ± %.2f%%\n', outName, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));
end



end
