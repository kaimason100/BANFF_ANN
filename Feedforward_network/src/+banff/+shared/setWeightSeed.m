function cleanup = setWeightSeed(seedValue)
% Set the global offset used by customWeights*_git. Default scripts use zero.
    setappdata(0, 'NAR_WEIGHT_SEED_BASE', seedValue);
    if nargout > 0, cleanup = onCleanup(@banff.shared.clearWeightSeed); end
end
