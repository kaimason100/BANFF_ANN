function trainTabular(task, seeds, cfg)
%TRAINTABULAR Shared bias-only training for tabular classification/regression.
% All numerical defaults and training-only normalisation follow the original
% release. Grouped splitting is the intentional change for duplicate records.
repoRoot = banff.root();
isRegression = any(strcmp(task,{'abalone','toyota'}));
if isRegression
    [X, targets, batchSize] = banff.shared.regressionDataset(repoRoot,task);
else
    [X, targets, batchSize] = banff.shared.classificationDataset(repoRoot,task);
end
for seed = seeds(:).'
    banff.shared.setWeightSeed(seed);
    cleanup = onCleanup(@banff.shared.clearWeightSeed); %#ok<NASGU>
    [trainIdx,valIdx,testIdx,splitPolicy] = banff.datasetSplit(X,targets,cfg.SplitSeed);
    featureStats = banff.shared.fitFeatureStandardization(X(:,trainIdx));
    banff.shared.assertFeatureStatsFromTraining(featureStats,X,trainIdx,task);
    XTrain = banff.shared.applyFeatureStandardization(X(:,trainIdx),featureStats);
    XVal = banff.shared.applyFeatureStandardization(X(:,valIdx),featureStats);
    banff.assertEffectivePartitions(banff.shared.applyFeatureStandardization(X,featureStats), ...
        {trainIdx,valIdx,testIdx});
    YTrain = targets(trainIdx); YVal = targets(valIdx);
    if isRegression
        targetStats.yMu = mean(YTrain);
        targetStats.yStd = std(YTrain);
        if targetStats.yStd == 0, targetStats.yStd = 1; end
        YTrain = (YTrain-targetStats.yMu)/targetStats.yStd;
        YVal = (YVal-targetStats.yMu)/targetStats.yStd;
        numOutputs = 1;
    else
        numOutputs = numel(categories(YTrain));
        assert(all(countcats(YTrain)>0),'A class is absent from the training partition.');
    end
    layers = [
        featureInputLayer(size(XTrain,1), 'Normalization','none')
        fullyConnectedLayer(cfg.Width,'WeightLearnRateFactor',0,'WeightsInitializer',@customWeights1_git,'BiasL2Factor',0,'BiasInitializer','zeros')
        tanhLayer
        fullyConnectedLayer(cfg.Width,'WeightLearnRateFactor',0,'WeightsInitializer',@customWeights2_git,'BiasL2Factor',0,'BiasInitializer','zeros')
        tanhLayer
        fullyConnectedLayer(numOutputs,'Name','output','WeightLearnRateFactor',0,'WeightsInitializer',@customWeights3_git,'BiasL2Factor',0,'BiasInitializer','zeros')];
    if isRegression, layers = [layers; regressionLayer];
    else, layers = [layers; softmaxLayer; classificationLayer]; end
    options = trainingOptions('adam', ...
        'MaxEpochs',cfg.MaxEpochs,'MiniBatchSize',batchSize, ...
        'InitialLearnRate',cfg.LearnRate,'Shuffle','every-epoch', ...
        'ValidationData',{XVal.',YVal},'ValidationFrequency',10, ...
        'OutputNetwork','best-validation-loss','Plots',banff.plotOption(cfg.Plots), ...
        'ExecutionEnvironment',cfg.ExecutionEnvironment,'Verbose',0, ...
        'OutputFcn',@(info) banff.progress(info,cfg.MaxEpochs,task));
    fprintf('[%s] training seed %d\n',task,seed);
    [net,info] = trainNetwork(XTrain.',YTrain,layers,options);
    split = struct('trainIdx',trainIdx,'valIdx',valIdx,'testIdx',testIdx,'policy',splitPolicy);
    metadata = struct('task',task,'seed',seed,'seedList',0:9,'width',cfg.Width, ...
        'normalization','train-fitted feature z-score','trainingOptions',cfg);
    metadata.datasetSHA256 = banff.fileSHA256(banff.datasetPath(task));
    modelPath = banff.networkPath(repoRoot,cfg.NetworkSet,task,seed);
    if isRegression
        metadata.normalization = 'train-fitted feature z-score with train-fitted target z-score';
        save(modelPath,'net','info','split','featureStats','targetStats','metadata','-v7.3');
    else
        save(modelPath,'net','info','split','featureStats','metadata','-v7.3');
    end
    clear cleanup
end
end
