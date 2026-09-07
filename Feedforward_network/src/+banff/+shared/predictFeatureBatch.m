function Y = predictFeatureBatch(net, X)
    assert(isnumeric(X) && ismatrix(X), 'Feature predictors must be a numeric 2-D features-by-observations array.');
    if isa(net, 'dlnetwork')
        net = dlupdate(@gather,net);
        Y = predict(net, dlarray(single(X), 'CB'));
    else
        Y = predict(net, single(X.'), 'ExecutionEnvironment','cpu');
    end
end
