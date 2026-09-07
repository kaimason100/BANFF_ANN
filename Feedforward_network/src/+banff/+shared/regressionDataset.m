function [X, y, batchSize] = regressionDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'abalone'
            S = load(fullfile(dataDir, 'abalone_dataset.mat')); D = S.data; D = replaceNonfiniteDatasetValues(D, task, 'loaded encoded matrix'); y = D(:, end); X = D(:, 1:end-1).'; batchSize = size(X, 2);
        case 'toyota'
            S = load(fullfile(dataDir, 'toyota_dataset.mat')); D = S.data; D = replaceNonfiniteDatasetValues(D, task, 'loaded encoded matrix'); y = D(:, 3); X = D(:, [1, 2, 4:size(D, 2)]).'; batchSize = min(4096, size(X, 2));
        otherwise
            error('Unknown regression task: %s', task);
    end
    X = replaceNonfiniteDatasetValues(X, 'current task', 'processed predictor or image values');
end
