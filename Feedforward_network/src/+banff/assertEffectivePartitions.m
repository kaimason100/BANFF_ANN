function assertEffectivePartitions(X, partitions)
%ASSERTEFFECTIVEPARTITIONS Check all pairs of feature-by-observation arrays.
for a=1:numel(partitions)
    for b=a+1:numel(partitions)
        banff.assertNoSharedInputs(single(X(:,partitions{a}).'), ...
            single(X(:,partitions{b}).'),sprintf('partitions %d/%d',a,b));
    end
end
end
