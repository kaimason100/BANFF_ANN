function clearWeightSeed()
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        rmappdata(0, 'NAR_WEIGHT_SEED_BASE');
    end
end
