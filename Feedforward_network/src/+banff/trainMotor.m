function trainMotor(weightSeed, cfg)
    weightCleanup = banff.shared.setWeightSeed(weightSeed); %#ok<NASGU>

repoRoot = banff.root();

episodeDir = fullfile(repoRoot,'outputs','training_data','motor_control');
if ~isfolder(episodeDir), mkdir(episodeDir); end
EPISODE_PATH = fullfile(episodeDir,sprintf('episodes_seed_%03d.mat',weightSeed));

rng(1, 'twister');

[A, B, K, dt, T, t, l1, l2] = banff.motor.lqrArmSetup();
numEpisodes = 250;
sampleInterval = 1;
inputsEpisodes = cell(numEpisodes, 1);
outputsEpisodes = cell(numEpisodes, 1);
targetEE = zeros(numEpisodes, 2);

for ep = 1:numEpisodes
    target = banff.motor.randomReachTarget(l1, l2);
    targetEE(ep, :) = target(:).';
    xTargetJoint = banff.motor.inverseKinematicsTarget(target, l1, l2);
    xSim = zeros(4, numel(t));
    epiInputs = [];
    epiOutputs = [];
    consecutiveCounter = 0;
    for k = 1:numel(t)-1
        eJoint = xSim(:, k) - xTargetJoint;
        u = -K * eJoint;
        xSim(:, k+1) = xSim(:, k) + dt*(A*xSim(:, k) + B*u);
        if mod(k, sampleInterval) == 0
            epiInputs = [epiInputs, banff.motor.encodeJointError(eJoint)]; %#ok<AGROW>
            epiOutputs = [epiOutputs, u]; %#ok<AGROW>
        end
        if banff.motor.endEffectorDistance(xSim(:, k+1), target, l1, l2) < 0.01
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
banff.motor.assertEpisodeSplit(trainEps, valEps, numEpisodes, 'LQR');
banff.motor.assertNoSharedTargets(targetEE(trainEps,:),targetEE(valEps,:),'LQR train/validation target leakage.');

XTrainAll = banff.motor.localConcat(inputsEpisodes, trainEps);
YTrainAll = banff.motor.localConcat(outputsEpisodes, trainEps);
XValAll = banff.motor.localConcat(inputsEpisodes, valEps);
YValAll = banff.motor.localConcat(outputsEpisodes, valEps);
assert(~isempty(XTrainAll) && ~isempty(XValAll), 'LQR train/validation samples are empty.');

mu_X = mean(XTrainAll, 2);
sigma_X = std(XTrainAll, 0, 2); sigma_X(sigma_X == 0) = 1;
mu_Y = mean(YTrainAll, 2);
sigma_Y = std(YTrainAll, 0, 2); sigma_Y(sigma_Y == 0) = 1;
banff.motor.assertStatsClose(mu_X, mean(XTrainAll, 2), 'LQR', 'input mean');
banff.motor.assertStatsClose(sigma_X, banff.motor.replaceZeroStd(std(XTrainAll, 0, 2)), 'LQR', 'input std');
banff.motor.assertStatsClose(mu_Y, mean(YTrainAll, 2), 'LQR', 'output mean');
banff.motor.assertStatsClose(sigma_Y, banff.motor.replaceZeroStd(std(YTrainAll, 0, 2)), 'LQR', 'output std');

numInputs = size(XTrainAll, 1);
XTrain = ((XTrainAll - mu_X) ./ (sigma_X * sqrt(numInputs))).';
YTrain = ((YTrainAll - mu_Y) ./ sigma_Y).';
XVal = ((XValAll - mu_X) ./ (sigma_X * sqrt(numInputs))).';
YVal = ((YValAll - mu_Y) ./ sigma_Y).';
banff.assertNoSharedInputs(single(XTrain),single(XVal),'Motor train/validation');
dataAudit = struct('XTrain',XTrainAll,'YTrain',YTrainAll,'XValidation',XValAll, ...
    'YValidation',YValAll,'targets',targetEE);

layers = [
    featureInputLayer(numInputs, 'Normalization', 'none', 'Name', 'input')
    fullyConnectedLayer(cfg.Width, 'Name', 'fc1', 'WeightsInitializer', @customWeights1_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')
    tanhLayer('Name', 'relu1')
    fullyConnectedLayer(cfg.Width, 'Name', 'fc2', 'WeightsInitializer', @customWeights2_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')
    tanhLayer('Name', 'relu2')
    fullyConnectedLayer(2, 'Name', 'fc3', 'WeightsInitializer', @customWeights3_git, 'WeightLearnRateFactor', 0, 'BiasInitializer', 'zeros')
    regressionLayer];

options = trainingOptions('adam', ...
    'Verbose', 0, 'MaxEpochs', cfg.MaxEpochs, 'InitialLearnRate', cfg.LearnRate, ...
    'L2Regularization', 0, 'Plots', banff.plotOption(cfg.Plots), 'ExecutionEnvironment', cfg.ExecutionEnvironment, 'MiniBatchSize', cfg.MiniBatchSize, ...
    'ValidationData', {XVal, YVal}, 'ValidationFrequency', 10, 'OutputNetwork', 'best-validation-loss', ...
    'OutputFcn', @(info) banff.progress(info,cfg.MaxEpochs,'Motor'));
[net, tr] = trainNetwork(XTrain, YTrain, layers, options);

numSims = 100;
testTargets = zeros(numSims, 2);
for sim = 1:numSims
    testTargets(sim, :) = banff.motor.randomReachTarget(l1, l2).';
end
banff.motor.assertNoSharedTargets(targetEE, testTargets, 'LQR train/test target leakage.');

testSequence.dt = dt; testSequence.T = T; testSequence.t = t; testSequence.l1 = l1; testSequence.l2 = l2;
testSequence.A = A; testSequence.B = B; testSequence.K = K;
testSequence.mu_X = mu_X; testSequence.sigma_X = sigma_X; testSequence.mu_Y = mu_Y; testSequence.sigma_Y = sigma_Y;
testSequence.numInputs = numInputs; testSequence.numSims = numSims; testSequence.sampleInterval = sampleInterval;
testSequence.trainTargets = targetEE(trainEps, :); testSequence.valTargets = targetEE(valEps, :); testSequence.testTargets = testTargets;
split.trainEps = trainEps; split.valEps = valEps;
metadata.task = 'LQR_two_link_arm'; metadata.normalization = 'train-episode-fitted input/output z-score';

metadata.seed = weightSeed; metadata.seedList = 0:9;
save(banff.networkPath(repoRoot,cfg.NetworkSet, 'LQR_two_link_arm', weightSeed), 'net', 'tr', 'mu_X', 'sigma_X', 'mu_Y', 'sigma_Y', 'split', 'testSequence', 'metadata', 'dataAudit', '-v7.3');
save(EPISODE_PATH, 'inputsEpisodes', 'outputsEpisodes', 'trainEps', 'valEps', 'targetEE', 'testSequence', '-v7.3');
end
