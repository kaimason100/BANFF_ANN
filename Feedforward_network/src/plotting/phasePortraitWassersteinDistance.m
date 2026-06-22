function distanceValue = phasePortraitWassersteinDistance(xPred, xTrue, options)
% phasePortraitWassersteinDistance  Sliced WD between 2-D phase portraits.
%
% The inputs are observations-by-state matrices. For systems with more than two
% states, the distance is averaged across all 2-D state-pair phase portraits.

if nargin < 3 || isempty(options)
    options = struct();
end
options = defaultPhaseWassersteinOptions(options);
xPred = double(xPred);
xTrue = double(xTrue);

if any(~isfinite(xPred(:))) || any(~isfinite(xTrue(:)))
    distanceValue = Inf;
    return
end

nStates = min(size(xPred, 2), size(xTrue, 2));
pairs = phasePortraitPairs(nStates);
pairScores = nan(size(pairs, 1), 1);
for p = 1:size(pairs, 1)
    predCloud = phasePortraitCloud(xPred(:, 1:nStates), pairs(p, :), options);
    trueCloud = phasePortraitCloud(xTrue(:, 1:nStates), pairs(p, :), options);
    if size(predCloud, 1) >= 2 && size(trueCloud, 1) >= 2
        pairScores(p) = slicedWasserstein2D(predCloud, trueCloud, options.NumProjections, options.TrimFraction);
    else
        pairScores(p) = Inf;
    end
end

finiteScores = pairScores(isfinite(pairScores));
if isempty(finiteScores)
    distanceValue = Inf;
else
    distanceValue = mean(finiteScores);
end
end

function options = defaultPhaseWassersteinOptions(options)
defaults = struct( ...
    'NumProjections', 128, ...
    'TrimFraction', 0.10, ...
    'Subsample', 5, ...
    'TransientFraction', 0.10, ...
    'MaxPoints', 1250);
names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(options, name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end
end

function pairs = phasePortraitPairs(nStates)
if nStates >= 2
    pairs = nchoosek(1:nStates, 2);
else
    pairs = [1 1];
end
end

function cloud = phasePortraitCloud(x, pair, options)
nRows = size(x, 1);
if nRows < 2
    cloud = zeros(0, 2);
    return
end

firstRow = 1 + floor(max(0, min(0.9, options.TransientFraction)) * nRows);
firstRow = min(max(firstRow, 1), nRows - 1);
x = x(firstRow:end, :);
x = x(1:max(1, round(options.Subsample)):end, :);

if pair(1) == pair(2)
    cloud = [x(1:end-1, pair(1)), x(2:end, pair(1))];
else
    cloud = x(:, pair);
end

cloud = cloud(all(isfinite(cloud), 2), :);
if size(cloud, 1) > options.MaxPoints
    idx = unique(round(linspace(1, size(cloud, 1), options.MaxPoints)));
    cloud = cloud(idx, :);
end
end

function d = slicedWasserstein2D(P, Q, numProjections, trimFraction)
if size(P, 1) < 2 || size(Q, 1) < 2
    d = Inf;
    return
end

theta = ((0:(numProjections-1)).' + 0.5) * pi / numProjections;
directions = [cos(theta), sin(theta)];
vals = nan(numProjections, 1);
for k = 1:numProjections
    p = sort(P * directions(k, :).');
    q = sort(Q * directions(k, :).');
    if numel(p) ~= numel(q)
        nInterp = max(numel(p), numel(q));
        t = linspace(0, 1, nInterp).';
        p = interp1(linspace(0, 1, numel(p)).', p, t, 'linear', 'extrap');
        q = interp1(linspace(0, 1, numel(q)).', q, t, 'linear', 'extrap');
    end
    vals(k) = mean((p - q).^2, 'omitnan');
end

vals = sort(vals(isfinite(vals)));
if isempty(vals)
    d = Inf;
    return
end
nTrim = floor(max(0, min(0.45, trimFraction)) * numel(vals));
if 2*nTrim < numel(vals)
    vals = vals((1+nTrim):(end-nTrim));
end
d = sqrt(mean(vals, 'omitnan'));
end
