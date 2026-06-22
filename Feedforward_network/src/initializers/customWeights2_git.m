function weights = customWeights2_git(sz)
%CUSTOMWEIGHTS2_GIT   Generate reproducible balanced ternary weight matrix.
%   WEIGHTS = CUSTOMWEIGHTS2_GIT(SZ) returns a matrix of size SZ
%   (a two-element vector [nRows, nCols]) whose entries are drawn from
%   {0, 1} or {0, –1} with equal column‐type probability and 2/3 chance
%   of the nonzero value. The result is scaled by 1/√16000.
%
%   Inputs:
%     sz   — 1×2 vector specifying [number of rows, number of columns].
%
%   Outputs:
%     weights — sz(1)×sz(2) matrix of balanced ternary weights.
%
%   Example:
%     W = customWeights2_git([512, 512]);
%
%   Notes:
%     • Full random pool has dimensions 16000×16000; sz must not exceed this.
%     • Columns are independently assigned type:
%         – Type "1" columns ∈ {0, 1}, with P(1)=2/3.
%         – Type "0" columns ∈ {0, –1}, with P(–1)=2/3.
%     • RNG seed fixed to 1 for deterministic output.

    % Input validation
    assert(isnumeric(sz) && numel(sz)==2 && all(sz>0) && all(mod(sz,1)==0), ...
        'sz must be a 1×2 vector of positive integers.');
    maxDim = 16e3;
    assert(all(sz <= maxDim), ...
        'sz dimensions must not exceed [%d, %d].', maxDim, maxDim);

    % Fix the random seed for reproducibility
    rng(9876 + getNarWeightSeedBase());

    % Define full pool size
    rows = maxDim;
    cols = maxDim;

    % Decide column types: 1 ⇒ {0,1}, 0 ⇒ {0,–1}
    col_types = randi([0, 1], 1, cols);

    % Preallocate full weight pool
    A = zeros(rows, cols);

    % Populate each column according to its type
    for c = 1:cols
        if col_types(c) == 1
            % Type {0,1}: element = 1 with probability 2/3, else 0
            A(:, c) = double(rand(rows, 1) < 2/3);
        else
            % Type {0,–1}: element = –1 with probability 2/3, else 0
            A(:, c) = -double(rand(rows, 1) < 2/3);
        end
    end

    % Scale by 1/sqrt(16000)
    W2 = A / sqrt(maxDim);

    % Return the requested submatrix
    weights = W2(1:sz(1), 1:sz(2));
end

function seedBase = getNarWeightSeedBase()
% Optional global seed offset for multi-seed release training.
    if isappdata(0, 'NAR_WEIGHT_SEED_BASE')
        seedBase = getappdata(0, 'NAR_WEIGHT_SEED_BASE');
    else
        seedBase = 0;
    end
end
