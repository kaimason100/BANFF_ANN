function assertDisjointSplits(trainIdx, valIdx, testIdx, n, label)
    trainIdx = trainIdx(:); valIdx = valIdx(:); testIdx = testIdx(:);
    assert(~isempty(trainIdx) && ~isempty(valIdx) && ~isempty(testIdx), '%s split is empty.', label);
    allIdx = [trainIdx; valIdx; testIdx];
    assert(all(allIdx >= 1 & allIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(isempty(intersect(trainIdx, testIdx)), '%s train/test leakage.', label);
    assert(isempty(intersect(valIdx, testIdx)), '%s validation/test leakage.', label);
    assert(numel(unique(allIdx)) == numel(allIdx), '%s duplicate split indices.', label);
end
