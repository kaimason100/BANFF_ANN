function assertPredictorPartitions(X, y, partitions)
%ASSERTPREDICTORPARTITIONS Reject index overlap, duplicate records and inputs
% shared between independent partitions. Different targets for the same input
% are permitted within one partition, as they are distinct observations.
y = y(:);
allIdx = [];
for k = 1:numel(partitions)
    idx = double(partitions{k}(:));
    assert(all(isfinite(idx) & idx == fix(idx) & idx >= 1 & idx <= size(X,2)), 'Invalid partition indices.');
    assert(isempty(intersect(allIdx, idx)), 'Partition indices overlap.');
    allIdx = [allIdx; idx]; %#ok<AGROW>
    rows = [double(X(:,idx).'), double(y(idx(:)))];
    assert(size(unique(rows, 'rows'),1) == numel(idx), 'Duplicate full records within a partition.');
    for j = 1:k-1
        assert(~any(ismember(X(:,idx).', X(:,partitions{j}).', 'rows')), ...
            'banff:DuplicateLeakage', 'Identical predictors cross saved partitions. Retrain with grouped splitting.');
    end
end
assert(numel(unique(allIdx)) == numel(allIdx), 'Duplicate partition indices.');
end
