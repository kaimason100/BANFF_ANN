function binIdx = safeColorBin(value, edges, nColours)
binIdx = discretize(value, edges);
if isnan(binIdx)
    if value <= edges(1)
        binIdx = 1;
    elseif value >= edges(end)
        binIdx = nColours;
    else
        binIdx = 1;
    end
end
binIdx = max(1, min(nColours, binIdx));
end
