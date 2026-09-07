function [trainIdx, valIdx, testIdx, policy] = datasetSplit(X, y, seed, hasTest)
%DATASETSPLIT Keep equal predictors in one partition; remove exact records.
% X is features by observations. Labels are used ONLY to identify duplicate
% records, never to assign partitions. Conflicting targets remain separate
% records in the same predictor group. The first exact record is retained.
% With no duplicates, preserve the original cvpartition calls and RNG state.
if nargin < 4, hasTest = true; end
assert(ismatrix(X) && isnumeric(X) && all(isfinite(X(:))), 'Predictors must be finite.');
assert(numel(y) == size(X, 2), 'Predictor and target counts differ.');
if iscategorical(y), assert(~any(isundefined(y)),'Undefined categorical targets.'); end
target = double(y(:));
assert(all(isfinite(target)), 'Targets must be finite.');
[~, keep] = unique([double(X.'), target], 'rows', 'stable');
keep = sort(keep);
[~, ~, group] = unique(X(:, keep).', 'rows', 'stable');
nGroups = max(group);
assert(nGroups >= 5, 'Too few distinct predictor groups to create held-out partitions.');
rng(double(seed), 'twister');
if hasTest
    cvTest = cvpartition(nGroups, 'HoldOut', 0.2);
    trainFullGroups = find(training(cvTest));
    testGroups = find(test(cvTest));
    cvVal = cvpartition(numel(trainFullGroups), 'HoldOut', 0.2);
    trainGroups = trainFullGroups(training(cvVal));
    valGroups = trainFullGroups(test(cvVal));
else
    cvVal = cvpartition(nGroups, 'HoldOut', 0.2);
    trainGroups = find(training(cvVal));
    valGroups = find(test(cvVal));
    testGroups = [];
end
trainIdx = keep(ismember(group, trainGroups));
valIdx = keep(ismember(group, valGroups));
testIdx = keep(ismember(group, testGroups));
policy = struct('name', 'exact-record-dedup-predictor-grouped-v1', ...
    'seed', double(seed), 'retainedRows', keep, ...
    'removedRows', setdiff((1:size(X, 2)).', keep), 'numPredictorGroups', nGroups);
banff.assertPredictorPartitions(X, y, {trainIdx, valIdx, testIdx});
if numel(keep) < size(X, 2) || nGroups < numel(keep)
    fprintf('Grouped split: removed %d exact records; %d distinct predictor groups across %d retained records.\n', ...
        size(X, 2)-numel(keep), nGroups, numel(keep));
end
end
