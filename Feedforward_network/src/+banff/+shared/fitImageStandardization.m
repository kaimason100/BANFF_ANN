function stats = fitImageStandardization(XTrain)
    XTrain = im2single(banff.shared.ensure4d(XTrain));
    stats.mu = mean(XTrain, 4);
    stats.sigma = std(XTrain, 0, 4);
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1 / sqrt(28*28);
end
