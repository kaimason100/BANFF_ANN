function idx = scoresToClassIndex(scores, nClasses, nObs)
    scores = double(scores);
    if size(scores, 1) == nObs && size(scores, 2) == nClasses
        [~, idx] = max(scores, [], 2);
    elseif size(scores, 1) == nClasses && size(scores, 2) == nObs
        [~, idx] = max(scores, [], 1);
        idx = idx(:);
    else
        error('Unexpected classification score size %s for %d classes and %d observations.', mat2str(size(scores)), nClasses, nObs);
    end
end
