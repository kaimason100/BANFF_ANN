function X = applyFeatureStandardization(X, stats)
    X = stats.scale * ((X - stats.mu) ./ stats.sigma);
    X = replaceNonfiniteDatasetValues(X, 'current task', 'processed predictor or image values');
end
