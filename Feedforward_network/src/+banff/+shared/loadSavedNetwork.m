function net = loadSavedNetwork(modelPath)
    S = load(modelPath);
    names = fieldnames(S);
    for candidate = ["net", "dlnet", "trainedNet", "netObj"]
        name = char(candidate);
        if isfield(S, name) && banff.shared.isNetworkLike(S.(name))
            net = S.(name);
            return
        end
    end
    for k = 1:numel(names)
        if banff.shared.isNetworkLike(S.(names{k}))
            net = S.(names{k});
            return
        end
    end
    error('No network-like variable found in %s.', modelPath);
end
