function X = applyImageStandardization(X, stats)
    X = im2single(banff.shared.ensure4d(X));
    X = stats.scale * ((X - stats.mu) ./ stats.sigma);
    X = replaceNonfiniteDatasetValues(X, 'current task', 'processed predictor or image values');
end
