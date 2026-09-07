function assertImageStatsFromTraining(stats, X, trainIdx, label)
    expected = banff.shared.fitImageStandardization(X(:, :, :, trainIdx));
    banff.shared.assertStatsClose(stats.mu, expected.mu, label, 'image mean');
    banff.shared.assertStatsClose(stats.sigma, expected.sigma, label, 'image standard deviation');
    banff.shared.assertStatsClose(stats.scale, expected.scale, label, 'image scale');
end
