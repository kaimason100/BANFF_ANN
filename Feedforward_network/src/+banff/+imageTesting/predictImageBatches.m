function Y = predictImageBatches(net, X, batchSize)
    X = banff.shared.ensure4d(X);
    nObs = size(X, 4);
    Y = [];
    for first = 1:batchSize:nObs
        last = min(first + batchSize - 1, nObs);
        XBatch = X(:, :, :, first:last);
        if banff.imageTesting.networkUsesImageInput(net)
            YBatch = banff.shared.numericOutput(predict(net, single(XBatch)));
        else
            XFlat = reshape(XBatch, [], numel(first:last));
            YBatch = banff.shared.numericOutput(banff.shared.predictFeatureBatch(net, XFlat));
        end
        if isempty(Y)
            if size(YBatch, 1) == numel(first:last)
                Y = zeros(nObs, size(YBatch, 2));
            else
                Y = zeros(size(YBatch, 1), nObs);
            end
        end
        if size(Y, 1) == nObs
            Y(first:last, :) = YBatch;
        else
            Y(:, first:last) = YBatch;
        end
    end
end
