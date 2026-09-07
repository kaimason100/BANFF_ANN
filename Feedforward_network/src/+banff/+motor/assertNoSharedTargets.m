function assertNoSharedTargets(trainTargets, testTargets, message)
    sharedTargets = intersect(round(trainTargets, 12), round(testTargets, 12), 'rows');
    assert(isempty(sharedTargets), message);
end
