function assertNoSharedInputs(A, B, label)
%ASSERTNOSHAREDINPUTS Reject shared effective predictors across partitions.
% Observations are rows. Call after the same scaling/cast used by the model.
assert(ismatrix(A) && ismatrix(B) && size(A,2)==size(B,2), ...
    'banff:InvalidPartition','%s: incompatible predictor dimensions.',label);
assert(~isempty(A) && ~isempty(B) && all(isfinite(A(:))) && all(isfinite(B(:))), ...
    'banff:InvalidPartition','%s: empty or non-finite predictors.',label);
assert(isempty(intersect(A,B,'rows')), 'banff:DuplicateLeakage', ...
    '%s: identical effective predictors occur in different partitions.',label);
end
