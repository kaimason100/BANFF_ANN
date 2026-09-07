function assertNoOverlap(trainIdx, valIdx, n, label)
    assert(~isempty(trainIdx) && ~isempty(valIdx), '%s train/validation split is empty.', label);
    assert(all(trainIdx >= 1 & trainIdx <= n) && all(valIdx >= 1 & valIdx <= n), '%s split indices out of range.', label);
    assert(isempty(intersect(trainIdx, valIdx)), '%s train/validation leakage.', label);
    assert(numel(unique([trainIdx(:); valIdx(:)])) == n, '%s duplicate or missing split indices.', label);
end
