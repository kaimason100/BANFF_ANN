function assertFeatureStatsFromTraining(stats, X, trainIdx, label)
    expected = banff.shared.fitFeatureStandardization(X(:, trainIdx));
    banff.shared.assertStatsClose(stats.mu, expected.mu, label, 'feature mean');
    banff.shared.assertStatsClose(stats.sigma, expected.sigma, label, 'feature standard deviation');
    banff.shared.assertStatsClose(stats.scale, expected.scale, label, 'feature scale');
end
