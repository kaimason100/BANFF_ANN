function [XTrain, YTrain, XTest, YTest, numOut] = mnistData(training, test, label)
    assert(isfield(training, 'images') && isfield(training, 'labels'), '%s missing training fields.', label);
    assert(isfield(test, 'images') && isfield(test, 'labels'), '%s missing test fields.', label);
    banff.shared.assertNoSharedImages(training.images, test.images, label);
    XTrain = banff.shared.ensure4d(training.images);
    XTest = banff.shared.ensure4d(test.images);
    YTrain = categorical(double(training.labels(:)) + 1);
    YTest = categorical(double(test.labels(:)) + 1);
    numOut = max(double(training.labels)) + 1;
    rows = reshape(XTest,[],size(XTest,4)).';
    [~,keep] = unique([double(rows),double(YTest(:))],'rows','stable');
    keep = sort(keep);
    XTest = XTest(:,:,:,keep); YTest = YTest(keep);
    assert(~any(isundefined(YTrain)) && ~any(isundefined(YTest)), 'Undefined image labels.');
end
