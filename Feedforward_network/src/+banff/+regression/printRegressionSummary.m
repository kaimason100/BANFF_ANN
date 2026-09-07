function printRegressionSummary(task, metrics)
    rmseMean = mean(metrics.rmse, 'omitnan');
    rmseSD = std(metrics.rmse, 0, 'omitnan');
    rMean = mean(metrics.r, 'omitnan');
    rSD = std(metrics.r, 0, 'omitnan');
    pMean = mean(metrics.p, 'omitnan');
    pSD = std(metrics.p, 0, 'omitnan');

    fprintf('%s seeded RMSE: mean %.4f | SD %.4f\n', task, rmseMean, rmseSD);
    fprintf('%s seeded Pearson r: mean %.4f | SD %.4f\n', task, rMean, rSD);
    fprintf('%s seeded Pearson p: mean %.4g | SD %.4g\n', task, pMean, pSD);
end
