function weights = customWeights1_git(sz)
%CUSTOMWEIGHTS1_GIT   Generate a reproducible weight matrix.
%   WEIGHTS = CUSTOMWEIGHTS1_GIT(SZ) returns a matrix of size SZ
%   (a two-element vector [nRows, nCols]) with entries drawn uniformly
%   from the interval [-2, 2]. The RNG is fixed for reproducibility.
%
%   Inputs:
%     sz   — 1×2 vector specifying [number of rows, number of columns].
%
%   Outputs:
%     weights — sz(1)×sz(2) matrix of uniformly sampled weights in [-2, 2].
%
%   Example:
%     W = customWeights1_git([512, 256]);
%
%   Notes:
%     • The full pool is 16 000×784; sz must not exceed these dimensions.
%     • RNG seed is set to 1 for deterministic behavior across runs.

    % Validate input
    assert(isnumeric(sz) && numel(sz)==2 && all(sz>0), ...
        'sz must be a 1×2 vector of positive integers.');
    assert(sz(1) <= 16000 && sz(2) <= 784, ...
        'sz dimensions must not exceed [16000, 784].');

    % Fix the random seed for reproducibility
    rng(1234 + getNarWeightSeedBase());

    % Build the full weight pool: 16000 rows × 784 columns, values ∼ U(–2,2)
    fullPool = 4 * rand(16000, 784) - 2;

    % Slice out the requested submatrix
    weights = fullPool(1:sz(1), 1:sz(2));
end

function seedBase = getNarWeightSeedBase()
% Optional global seed offset for multi-seed release training.
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        seedBase = getappdata(0, 'NAR_WEIGHT_SEED_BASE');
    else
        seedBase = 0;
    end
end
