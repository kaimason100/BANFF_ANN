function testClassification(taskName, SEEDS)
%% Aggregate seeded classification test results
repoRoot = banff.root();


tasks = string(taskName);
for task = tasks
    task = char(task);
    [X, labels, ~] = banff.shared.classificationDataset(repoRoot, task);
    [expectedTrainIdx, expectedValIdx, expectedTestIdx] = banff.datasetSplit(X, labels, 42);
    acc = nan(numel(SEEDS), 1);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded_grouped', task, sprintf('%s_seed_%03d_network.mat', task, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assert(isfield(S.metadata,'datasetSHA256') && strcmp(S.metadata.datasetSHA256, ...
            banff.fileSHA256(banff.datasetPath(task))), 'Dataset provenance mismatch; retrain with this branch.');
        assert(isfield(S, 'split'), 'Seeded model lacks split metadata: %s', modelPath);
        banff.shared.assertSeedMetadata(S, task, seedValue, SEEDS);
        assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
        banff.shared.assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), task);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
            '%s seed %d has a different split. Retrain with this branch; legacy splits are not accepted.', task, seedValue);
        banff.shared.assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, task);
        banff.assertPredictorPartitions(X, labels, {S.split.trainIdx,S.split.valIdx,S.split.testIdx});
        banff.assertEffectivePartitions(banff.shared.applyFeatureStandardization(X,S.featureStats), ...
            {S.split.trainIdx,S.split.valIdx,S.split.testIdx});
        net = banff.shared.loadSavedNetwork(modelPath);
        XTest = banff.shared.applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats); YTest = labels(S.split.testIdx);
        YTrain = labels(S.split.trainIdx);
        scores = banff.shared.numericOutput(banff.shared.predictFeatureBatch(net, XTest));
        predIdx = banff.classification.scoresToClassIndex(scores, numel(categories(YTrain)), numel(YTest));
        yHat = categorical(predIdx(:), 1:numel(categories(YTrain)), categories(YTrain));
        acc(i) = mean(yHat(:) == YTest(:));
        fprintf('%s seed %d accuracy: %.2f%%\n', task, seedValue, 100*acc(i));
    end
    fprintf('%s seeded accuracy mean ± SD: %.2f ± %.2f%%\n', task, 100*mean(acc, 'omitnan'), 100*std(acc, 0, 'omitnan'));
end



end
