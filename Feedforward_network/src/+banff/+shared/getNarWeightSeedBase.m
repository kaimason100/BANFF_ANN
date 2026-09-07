function seedBase = getNarWeightSeedBase()
% Optional global seed offset for multi-seed release training.
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        seedBase = getappdata(0, 'NAR_WEIGHT_SEED_BASE');
    else
        seedBase = 0;
    end
end
