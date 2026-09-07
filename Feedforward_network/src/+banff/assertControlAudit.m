function assertControlAudit(S, task)
%ASSERTCONTROLAUDIT Recheck saved control partitions and training-only scaling.
assert(isfield(S,'dataAudit'),'banff:MissingAudit', ...
    '%s lacks partition audit data. Retrain with this branch.',task);
A = S.dataAudit;
if strcmp(task,'Pong')
    banff.assertPredictorPartitions(A.features.',A.labels,{S.split.trainIdx,S.split.valIdx});
    banff.assertNoSharedInputs(single(A.features(S.split.trainIdx,:)), ...
        single(A.features(S.split.valIdx,:)),'Pong train/validation');
else
    banff.motor.assertEpisodeSplit(S.split.trainEps,S.split.valEps,size(A.targets,1),task);
    banff.motor.assertNoSharedTargets(A.targets(S.split.trainEps,:),A.targets(S.split.valEps,:), ...
        'Motor train/validation targets overlap.');
    banff.shared.assertStatsClose(S.mu_X,mean(A.XTrain,2),task,'input mean');
    banff.shared.assertStatsClose(S.sigma_X,banff.motor.replaceZeroStd(std(A.XTrain,0,2)),task,'input SD');
    banff.shared.assertStatsClose(S.mu_Y,mean(A.YTrain,2),task,'output mean');
    banff.shared.assertStatsClose(S.sigma_Y,banff.motor.replaceZeroStd(std(A.YTrain,0,2)),task,'output SD');
    X = (A.XTrain-S.mu_X)./(S.sigma_X*sqrt(size(A.XTrain,1)));
    V = (A.XValidation-S.mu_X)./(S.sigma_X*sqrt(size(A.XTrain,1)));
    banff.assertNoSharedInputs(single(X.'),single(V.'),'Motor train/validation');
    assert(isequal(S.testSequence.trainTargets,A.targets(S.split.trainEps,:)) && ...
        isequal(S.testSequence.valTargets,A.targets(S.split.valEps,:)), ...
        'banff:SourceMismatch','Motor saved target provenance differs.');
end
end
