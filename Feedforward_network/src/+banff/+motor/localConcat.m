function M = localConcat(C, idx)
    M = [];
    for ii = 1:numel(idx)
        ci = C{idx(ii)};
        if ~isempty(ci), M = [M, ci]; end %#ok<AGROW>
    end
end
