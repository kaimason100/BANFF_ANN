function modelPath = seededModelPath(repoRoot, taskName, seedValue)
    folderName = matlab.lang.makeValidName(char(taskName));
    modelPath = fullfile(repoRoot, 'trained_networks', 'seeded_grouped', folderName, sprintf('%s_seed_%03d_network.mat', folderName, double(seedValue)));
    assert(isfile(modelPath), 'Missing seeded model: %s', modelPath);
end
