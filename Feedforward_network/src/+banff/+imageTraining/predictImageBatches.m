function Y = predictImageBatches(net, X, batchSize)
    nObs = size(X, 4);
    Y = [];
    for first = 1:batchSize:nObs
        last = min(first + batchSize - 1, nObs);
        YBatch = banff.shared.numericOutput(predict(net, X(:, :, :, first:last)));
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
