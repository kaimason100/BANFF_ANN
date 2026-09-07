function [X, labels, batchSize] = classificationDataset(repoRoot, task)
    dataDir = fullfile(repoRoot, 'data');
    switch task
        case 'iris'
            load fisheriris
            labels = categorical(species);
            X = meas.';
            batchSize = 150;
        case 'breast_cancer'
            S = load(fullfile(dataDir, 'breast_cancer_dataset.mat')); [labels, X] = banff.shared.breastCancerData(S.data); batchSize = 256;
        case 'car_quality'
            S = load(fullfile(dataDir, 'car_dataset.mat')); D = S.data; labels = categorical(D(:, 7)); D = replaceNonfiniteDatasetValues(D, task, 'loaded encoded matrix'); X = D(:, 1:6).'; batchSize = 256;
        case 'mushroom'
            S = load(fullfile(dataDir, 'mushroom_dataset.mat')); D = S.data; labels = categorical(D(:, 1)); D = replaceNonfiniteDatasetValues(D, task, 'loaded encoded matrix'); X = D(:, 2:end).'; batchSize = 4096;
        otherwise
            error('Unknown classification task: %s', task);
    end
    X = replaceNonfiniteDatasetValues(X, 'current task', 'processed predictor or image values');
end
