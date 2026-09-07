function tf = isNetworkLike(v)
    c = lower(class(v));
    tf = contains(c, 'network') || strcmp(c, 'dlnetwork') || strcmp(c, 'seriesnetwork') || strcmp(c, 'dagnetwork');
end
