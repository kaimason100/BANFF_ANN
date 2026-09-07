function addReleasePaths(repoRoot)
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(genpath(fullfile(repoRoot, 'examples')));
    addpath(genpath(fullfile(repoRoot, 'src', 'dynamical_systems')));
end
