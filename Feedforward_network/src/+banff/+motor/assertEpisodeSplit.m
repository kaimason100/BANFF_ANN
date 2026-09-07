function assertEpisodeSplit(trainEps, valEps, n, label)
    assert(~isempty(trainEps) && ~isempty(valEps), '%s episode split is empty.', label);
    assert(all(trainEps >= 1 & trainEps <= n) && all(valEps >= 1 & valEps <= n), '%s episode indices out of range.', label);
    assert(isempty(intersect(trainEps, valEps)), '%s train/validation episode leakage.', label);
    assert(numel(unique([trainEps(:); valEps(:)])) == n, '%s duplicate or missing episode indices.', label);
end
