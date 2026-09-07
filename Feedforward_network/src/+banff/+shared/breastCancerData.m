function [labels, X] = breastCancerData(data)
    if istable(data)
        labels = categorical(string(data{:, 2})); X = data{:, 3:end}.';
    elseif iscell(data)
        labels = categorical(string(data(:, 2))); X = cellfun(@double, data(:, 3:end)).';
    elseif isnumeric(data)
        [~, ~, g] = unique(data(:, 2)); labels = categorical(g); X = data(:, 3:end).';
    elseif isstruct(data) && isfield(data, 'features') && isfield(data, 'labels')
        X = data.features.'; labels = categorical(data.labels(:));
    else
        error('Unsupported breast-cancer data format.');
    end
    X = replaceNonfiniteDatasetValues(X, 'current task', 'processed predictor or image values');
end
