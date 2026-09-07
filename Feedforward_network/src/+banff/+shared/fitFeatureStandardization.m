function stats = fitFeatureStandardization(XTrain)
    stats.mu = mean(XTrain, 2);
    stats.sigma = std(XTrain, 0, 2);
    stats.sigma(~isfinite(stats.sigma) | stats.sigma == 0) = 1;
    stats.scale = 1 / sqrt(size(XTrain, 1));
end
