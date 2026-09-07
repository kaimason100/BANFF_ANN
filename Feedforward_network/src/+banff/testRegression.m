function testRegression(taskName, SEEDS)
%% Aggregate seeded regression test results
repoRoot = banff.root();


tasks = string(taskName);
for task = tasks
    task = char(task);
    [X, y, ~] = banff.shared.regressionDataset(repoRoot, task);
    [expectedTrainIdx, expectedValIdx, expectedTestIdx] = banff.datasetSplit(X, y, banff.splitSeed(task));
    metrics = struct('rmse', [], 'r', [], 'p', []);
    for i = 1:numel(SEEDS)
        seedValue = SEEDS(i);
        modelPath = fullfile(repoRoot, 'trained_networks', 'seeded_grouped', task, sprintf('%s_seed_%03d_network.mat', task, seedValue));
        assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
        S = load(modelPath);
        assert(isfield(S.metadata,'datasetSHA256') && strcmp(S.metadata.datasetSHA256, ...
            banff.fileSHA256(banff.datasetPath(task))), 'Dataset provenance mismatch; retrain with this branch.');
        banff.shared.assertSeedMetadata(S, task, seedValue, SEEDS);
        assert(isfield(S, 'featureStats'), 'Seeded model lacks train-fitted feature normalization metadata: %s', modelPath);
        assert(isfield(S, 'targetStats') && isfield(S.targetStats, 'yMu') && isfield(S.targetStats, 'yStd'), 'Seeded regression model lacks target normalization metadata: %s', modelPath);
        banff.shared.assertDisjointSplits(S.split.trainIdx, S.split.valIdx, S.split.testIdx, size(X, 2), task);
        assert(isequal(S.split.trainIdx(:), expectedTrainIdx(:)) && isequal(S.split.valIdx(:), expectedValIdx(:)) && isequal(S.split.testIdx(:), expectedTestIdx(:)), ...
            '%s seed %d has a different split. Retrain with this branch; legacy splits are not accepted.', task, seedValue);
        banff.shared.assertFeatureStatsFromTraining(S.featureStats, X, S.split.trainIdx, task);
        banff.assertPredictorPartitions(X, y, {S.split.trainIdx,S.split.valIdx,S.split.testIdx});
        banff.assertEffectivePartitions(banff.shared.applyFeatureStandardization(X,S.featureStats), ...
            {S.split.trainIdx,S.split.valIdx,S.split.testIdx});
        net = banff.shared.loadSavedNetwork(modelPath);
        XTest = banff.shared.applyFeatureStandardization(X(:, S.split.testIdx), S.featureStats);
        YpredN = banff.shared.numericOutput(banff.shared.predictFeatureBatch(net, XTest));
        Ypred = YpredN(:) * S.targetStats.yStd + S.targetStats.yMu;
        Ytrue = double(y(S.split.testIdx));
        assert(numel(Ypred)==numel(Ytrue) && all(isfinite(Ypred)), 'Invalid regression predictions.');
        expectedStd = std(y(S.split.trainIdx)); if expectedStd == 0, expectedStd = 1; end
        banff.shared.assertStatsClose(S.targetStats.yMu,mean(y(S.split.trainIdx)),task,'target mean');
        banff.shared.assertStatsClose(S.targetStats.yStd,expectedStd,task,'target standard deviation');
        err = Ypred(:) - Ytrue(:);
        metrics.rmse(i,1) = sqrt(mean(err.^2, 'omitnan'));
        [metrics.r(i,1), metrics.p(i,1)] = banff.regression.corrStats(Ypred(:), Ytrue(:));
    end
    banff.regression.printRegressionSummary(task, metrics);
end



end
