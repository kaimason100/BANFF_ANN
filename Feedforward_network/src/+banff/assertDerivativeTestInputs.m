function assertDerivativeTestInputs(S, x0)
%ASSERTDERIVATIVETESTINPUTS Reconstruct fixed partitions without changing RNG.
M = S.metadata;
keys = {'sampleRandomSeed','validationRandomSeed','randomSamplesPerEpoch','validationSampleCount','stateRandomRange'};
for k=1:numel(keys)
    assert(isfield(M,keys{k}),'banff:MissingAudit','Missing derivative-field sample provenance: %s.',keys{k});
end
lo=double(M.stateRandomRange(1)); hi=double(M.stateRandomRange(2));
d=numel(S.stateStats.mu);
trainStream=RandStream('mt19937ar','Seed',double(M.sampleRandomSeed));
valStream=RandStream('mt19937ar','Seed',double(M.validationRandomSeed));
X=lo+(hi-lo)*rand(trainStream,double(M.randomSamplesPerEpoch),d);
V=lo+(hi-lo)*rand(valStream,double(M.validationSampleCount),d);
banff.assertNoSharedInputs(single(X),single(V),'Derivative train/validation');
banff.assertNoSharedInputs(single([X;V]),single(x0(:).'),'Derivative development/test IC');
if isfield(S,'dataAudit')
    assert(isequal(X,S.dataAudit.XTrain) && isequal(V,S.dataAudit.XValidation), ...
        'banff:SourceMismatch','Saved derivative samples disagree with their provenance.');
end
end
