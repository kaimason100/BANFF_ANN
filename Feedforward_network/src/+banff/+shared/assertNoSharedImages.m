function assertNoSharedImages(trainImages, testImages, label)
    trainImages = banff.shared.ensure4d(trainImages);
    testImages = banff.shared.ensure4d(testImages);
    trainRows = reshape(trainImages, [], size(trainImages, 4)).';
    testRows = reshape(testImages, [], size(testImages, 4)).';
    assert(isempty(intersect(trainRows, testRows, 'rows')), '%s exact duplicate images across train/test.', label);
end
