function trainPong(weightSeed, cfg)
    weightCleanup = banff.shared.setWeightSeed(weightSeed); %#ok<NASGU>

repoRoot = banff.root();


rng(1, 'twister');
ballSpeed = 0.02;
maxBounceAngle = pi/4;
numSamples = 1000;
[features, labels] = banff.pong.generateTrainingData(numSamples, ballSpeed);

numTrain = round(0.8 * numSamples);
trainIdx = 1:numTrain;
valIdx = (numTrain+1):numSamples;
banff.pong.assertNoOverlap(trainIdx, valIdx, numSamples, 'Pong');
banff.assertPredictorPartitions(features.',labels,{trainIdx,valIdx});
trainFeatures = features(trainIdx, :); trainLabels = labels(trainIdx);
valFeatures = features(valIdx, :); valLabels = labels(valIdx);
banff.assertNoSharedInputs(single(trainFeatures),single(valFeatures),'Pong train/validation');
dataAudit = struct('features',features,'labels',labels);

layers = banff.pong.createPongNetwork(cfg.Width);

options = trainingOptions('adam', ...
    'MaxEpochs', cfg.MaxEpochs, 'MiniBatchSize', cfg.MiniBatchSize, 'ValidationData', {valFeatures, valLabels}, ...
    'ValidationFrequency', 10, 'OutputNetwork', 'best-validation-loss', ...
    'OutputFcn', @(info) banff.progress(info,cfg.MaxEpochs,'Pong'), ...
    'Plots', banff.plotOption(cfg.Plots), 'ExecutionEnvironment', cfg.ExecutionEnvironment, 'Verbose', 0, 'InitialLearnRate', cfg.LearnRate);
[net, tr] = trainNetwork(trainFeatures, trainLabels, layers, options);

scores = banff.shared.numericOutput(predict(net, valFeatures));
predictedIdx = banff.pong.scoresToClassIndex(scores, 3, numel(valLabels));
accuracy = mean(predictedIdx(:) == grp2idx(valLabels(:)));
fprintf('Pong validation accuracy: %.2f%%\n', 100*accuracy);

split.trainIdx = trainIdx; split.valIdx = valIdx;
metadata.task = 'Pong'; metadata.normalization = 'none'; metadata.ballSpeed = ballSpeed; metadata.maxBounceAngle = maxBounceAngle;
metadata.seed = weightSeed; metadata.seedList = 0:9;
save(banff.networkPath(repoRoot,cfg.NetworkSet, 'Pong', weightSeed), 'net', 'tr', 'split', 'metadata', 'dataAudit', 'ballSpeed', 'maxBounceAngle', '-v7.3');
end
