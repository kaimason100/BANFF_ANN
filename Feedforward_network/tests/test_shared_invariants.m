function tests = test_shared_invariants
%TEST_SHARED_INVARIANTS Fast MATLAB tests; no network training or GPU required.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'setup.m'));
testCase.TestData.root = root;
end

function testOriginalSplitAndRNGWithoutDuplicates(testCase)
X = reshape(1:600,3,200); y = (1:200).';
[a,b,c] = banff.datasetSplit(X,y,42);
actualRNG = rng;
rng(42,'twister');
cv = cvpartition(200,'HoldOut',0.2);
fullTrain = find(training(cv)); expectedTest = find(test(cv));
cvVal = cvpartition(numel(fullTrain),'HoldOut',0.2);
verifyEqual(testCase,a,fullTrain(training(cvVal)));
verifyEqual(testCase,b,fullTrain(test(cvVal)));
verifyEqual(testCase,c,expectedTest);
verifyEqual(testCase,actualRNG,rng);
end

function testDuplicatesAndConflictingTargets(testCase)
X = [1:100,1,1,2]; y = [(1:100).';1;999;2];
[a,b,c,p] = banff.datasetSplit(X,y,42);
verifyEqual(testCase,sort(p.removedRows),[101;103]);
parts = {a,b,c};
for k = 1:3
    verifyEqual(testCase,ismember(1,parts{k}),ismember(102,parts{k}));
end
banff.assertPredictorPartitions(X,y,parts);
verifyEqual(testCase,numel([a;b;c]),101);
end

function testImageSplitPreservesOrderAndRNG(testCase)
images = reshape(uint8(1:100),1,1,1,100);
labels = categorical(mod((1:100).',10));
[a,b] = banff.imageSplit(images,labels,42);
actualRNG = rng;
rng(42,'twister'); order = randperm(100);
verifyEqual(testCase,a,order(1:80));
verifyEqual(testCase,b,order(81:end));
verifyEqual(testCase,actualRNG,rng);
end

function testPublicationSchedule(testCase)
S = banff.publicationSchedule();
verifyEqual(testCase,height(S),250);
verifyEqual(testCase,numel(unique(S.Task)),24);
initial = S.Stage == "initial";
verifyEqual(testCase,sum(initial),240);
verifyEqual(testCase,sum(S.DerivativeField & initial),90);
exceptions = S(S.LearnRate==0.005,:);
verifyEqual(testCase,sort(exceptions.Task+"/"+string(exceptions.Seed)), ...
    sort(["MO5/3";"MO7/8";"MO13/3";"Rikitake/0";"Rikitake/4"]));
continuation = S(S.Stage=="continuation",:);
verifyEqual(testCase,continuation.Task,repmat("MO0",10,1));
verifyEqual(testCase,continuation.Epochs,repmat(50000,10,1));
verifyEqual(testCase,continuation.LearnRate,repmat(0.001,10,1));
end

function testInputLeakageIsRejected(testCase)
verifyError(testCase,@() banff.assertPredictorPartitions([1,2,1],[0;0;0],{[1,2],3}), ...
    'banff:DuplicateLeakage');
end

function testEffectivePrecisionCollisionIsRejected(testCase)
verifyError(testCase,@() banff.assertNoSharedInputs(single(1),single(1+1e-9),'precision'), ...
    'banff:DuplicateLeakage');
end

function testDerivativeIdentity(testCase)
S.metadata = struct('task','Lorenz','seed',2, ...
    'trainingSamplingMode','state-random-derivative-field');
banff.assertDerivativeModelIdentity(S,'Lorenz',2);
verifyError(testCase,@() banff.assertDerivativeModelIdentity(S,'MO0',2),'banff:ModelIdentity');
verifyError(testCase,@() banff.assertDerivativeModelIdentity(S,'Lorenz',3),'banff:ModelIdentity');
S.metadata.sourceNetworkSet = 'other';
base = fullfile(tempdir,'base','Lorenz','Lorenz_seed_002_network.mat');
verifyError(testCase,@() banff.assertDerivativeModelIdentity(S,'Lorenz',2,base,true),'banff:ModelIdentity');
S.metadata.sourceNetworkSet = 'base';
verifyWarning(testCase,@() banff.assertDerivativeModelIdentity(S,'Lorenz',2,base,true), ...
    'banff:UnverifiedContinuationSource');
S.metadata.sourceNetworkSHA256 = 'deliberately-invalid';
verifyError(testCase,@() banff.assertDerivativeModelIdentity(S,'Lorenz',2,base,true),'banff:ModelIdentity');
end

function testTabularBatchOptions(testCase)
for task = {'iris','breast_cancer','car_quality','mushroom','abalone','toyota'}
    cfg = banff.settings(task{1});
    verifyEmpty(testCase,cfg.MiniBatchSize);
    cfg = banff.settings(task{1},struct('MiniBatchSize',23));
    verifyEqual(testCase,cfg.MiniBatchSize,23);
end
end

function testPublicationWDUsesMetricGrid(testCase)
t = (0:0.01:2).';
X = [sin(t),cos(t),sin(2*t)];
D.phaseWD = struct('time',t,'true',X,'output',X+0.1,'outputDt',0.01);
% Deliberately different plotting data must not affect the reported metric.
D.true = zeros(3); D.output = ones(3);
options = struct('NumProjections',128,'TrimFraction',0.1, ...
    'Subsample',5,'TransientFraction',0.1,'MaxPoints',1250);
verifyEqual(testCase,publicationPhaseWassersteinDistance(D,options), ...
    phasePortraitWassersteinDistance(X+0.1,X,options));
D.phaseWD.time(2) = 0.009;
verifyError(testCase,@() publicationPhaseWassersteinDistance(D,options),'banff:InvalidTestGrid');
verifyError(testCase,@() publicationPhaseWassersteinDistance(struct(),options),'banff:MissingTestGrid');
end
