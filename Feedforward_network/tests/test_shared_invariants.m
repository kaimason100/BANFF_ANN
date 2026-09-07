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
