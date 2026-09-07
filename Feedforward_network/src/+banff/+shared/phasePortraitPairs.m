function pairs = phasePortraitPairs(nStates)
if nStates >= 2
    pairs = nchoosek(1:nStates, 2);
else
    pairs = [1 1];
end
end
