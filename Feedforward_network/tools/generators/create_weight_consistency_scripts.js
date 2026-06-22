const fs = require("fs");
const path = require("path");
const os = require("os");
const childProcess = require("child_process");

const repo = process.cwd();
const template = path.join(repo, "examples", "classification", "Iris_classification.mlx");

function esc(text) {
  return text.replaceAll("]]>", "]]]]><![CDATA[>");
}

function makeMlx(outPath, code) {
  const work = fs.mkdtempSync(path.join(os.tmpdir(), "mlx-weights-"));
  childProcess.execFileSync("unzip", ["-q", template, "-d", work]);
  fs.writeFileSync(
    path.join(work, "matlab", "document.xml"),
    `<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body><w:p><w:pPr><w:pStyle w:val="code"/></w:pPr><w:r><w:t><![CDATA[${esc(code)}]]></w:t></w:r></w:p></w:body></w:document>`
  );
  fs.writeFileSync(
    path.join(work, "matlab", "output.xml"),
    `<?xml version="1.0" encoding="UTF-8"?><embeddedOutputs><metaData><evaluationState>manual</evaluationState><layoutState>code</layoutState><outputStatus>ready</outputStatus></metaData><outputArray type="array"/><regionArray type="array"/></embeddedOutputs>`
  );
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.rmSync(outPath, { force: true });
  childProcess.execFileSync("zip", ["-qr", outPath, "."], { cwd: work });
  fs.rmSync(work, { recursive: true, force: true });
}

const helpers = String.raw`
function repoRoot = locateRepoRoot()
% Return the release root when this live script is run from any folder.
    starts = string.empty;
    scriptPath = mfilename('fullpath');
    if ~isempty(scriptPath)
        starts(end+1) = string(fileparts(scriptPath)); %#ok<AGROW>
    end
    stack = dbstack('-completenames');
    for k = 1:numel(stack)
        if isfield(stack(k), 'file') && ~isempty(stack(k).file)
            starts(end+1) = string(fileparts(stack(k).file)); %#ok<AGROW>
        end
    end
    starts(end+1) = string(pwd);
    for s = starts
        candidate = char(s);
        while ~isempty(candidate)
            if isfolder(fullfile(candidate, 'trained_networks')) && isfolder(fullfile(candidate, 'examples'))
                repoRoot = candidate;
                return
            end
            parent = fileparts(candidate);
            if strcmp(parent, candidate), break; end
            candidate = parent;
        end
    end
    error('Could not locate release root. Run from inside the release folder.');
end

function addReleasePaths(repoRoot)
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(genpath(fullfile(repoRoot, 'examples')));
    addpath(genpath(fullfile(repoRoot, 'src', 'dynamical_systems')));
end

function records = loadWeightRecords(files)
% Load each saved network and extract encoder, hidden, and decoder weights.
    records = struct('label', {}, 'file', {}, 'weights', {});
    for i = 1:numel(files)
        modelPath = fullfile(files(i).folder, files(i).name);
        [net, label] = loadSavedNetwork(modelPath);
        weights = extractWeightTriplet(net, label);
        records(end+1).label = label; %#ok<AGROW>
        records(end).file = modelPath;
        records(end).weights = weights;
        fprintf('Loaded %-40s encoder %s | hidden %s | decoder %s\n', ...
            label, mat2str(size(weights.encoder)), mat2str(size(weights.hidden)), mat2str(size(weights.decoder)));
    end
end

function [net, label] = loadSavedNetwork(modelPath)
    assert(isfile(modelPath), 'Missing saved network: %s', modelPath);
    S = load(modelPath);
    label = erase(string(getFileName(modelPath)), ".mat");
    for candidate = ["net", "dlnet", "trainedNet", "netObj"]
        name = char(candidate);
        if isfield(S, name) && isNetworkLike(S.(name))
            net = S.(name);
            return
        end
    end
    names = fieldnames(S);
    for k = 1:numel(names)
        if isNetworkLike(S.(names{k}))
            net = S.(names{k});
            return
        end
    end
    error('%s has no loadable network object.', modelPath);
end

function name = getFileName(filePath)
    [~, base, ext] = fileparts(filePath);
    name = [base ext];
end

function tf = isNetworkLike(v)
    c = lower(class(v));
    tf = contains(c, 'network') || strcmp(c, 'dlnetwork') || strcmp(c, 'seriesnetwork') || strcmp(c, 'dagnetwork');
end

function weights = extractWeightTriplet(net, label)
% The release networks are expected to have encoder, hidden-to-hidden, decoder FC weights.
    fcWeights = extractFullyConnectedWeights(net);
    assert(numel(fcWeights) >= 3, '%s has %d fully connected weight matrices; expected at least 3.', label, numel(fcWeights));
    weights.encoder = fcWeights{1};
    weights.hidden = fcWeights{2};
    weights.decoder = fcWeights{end};
end

function fcWeights = extractFullyConnectedWeights(net)
    fcWeights = {};
    if isprop(net, 'Learnables')
        T = net.Learnables;
        params = string(T.Parameter);
        rows = find(params == "Weights");
        for i = 1:numel(rows)
            fcWeights{end+1} = numericMatrix(T.Value{rows(i)}); %#ok<AGROW>
        end
        if ~isempty(fcWeights), return; end
    end
    assert(isprop(net, 'Layers'), 'Network object has neither Learnables nor Layers.');
    layers = net.Layers;
    for i = 1:numel(layers)
        if isprop(layers(i), 'Weights') && ~isempty(layers(i).Weights)
            fcWeights{end+1} = numericMatrix(layers(i).Weights); %#ok<AGROW>
        end
    end
end

function A = numericMatrix(A)
    if isa(A, 'dlarray'), A = extractdata(A); end
    try, A = gather(A); catch, end
    A = double(A);
    assert(ismatrix(A), 'Expected a matrix weight tensor, got size %s.', mat2str(size(A)));
end

function compareAllRecords(records, titleText)
    assert(numel(records) >= 2, '%s needs at least two networks to compare.', titleText);
    tol = 0;
    fprintf('\n%s\n', titleText);
    for i = 1:numel(records)-1
        for j = i+1:numel(records)
            comparePair(records(i), records(j), tol);
        end
    end
    fprintf('%s passed for %d networks.\n', titleText, numel(records));
end

function comparePair(a, b, tol)
    encDiff = compareCommonMatrix(a.weights.encoder, b.weights.encoder, 'encoder', a.label, b.label, tol);
    hiddenDiff = compareExactMatrix(a.weights.hidden, b.weights.hidden, 'hidden-to-hidden', a.label, b.label, tol);
    decDiff = compareCommonMatrix(a.weights.decoder, b.weights.decoder, 'decoder', a.label, b.label, tol);
    fprintf('%s vs %s | encoder %.3g | hidden %.3g | decoder %.3g\n', a.label, b.label, encDiff, hiddenDiff, decDiff);
end

function maxDiff = compareCommonMatrix(A, B, partName, labelA, labelB, tol)
% Compare the shared leading submatrix when input/output counts differ.
    nRows = min(size(A, 1), size(B, 1));
    nCols = min(size(A, 2), size(B, 2));
    assert(nRows > 0 && nCols > 0, '%s comparison between %s and %s has an empty common submatrix.', partName, labelA, labelB);
    maxDiff = maxAbsDiffLeading(A, B, nRows, nCols);
    assert(maxDiff <= tol, '%s weights differ between %s and %s over the common %d-by-%d submatrix. Max abs diff %.17g.', ...
        partName, labelA, labelB, nRows, nCols, maxDiff);
end

function maxDiff = compareExactMatrix(A, B, partName, labelA, labelB, tol)
% Hidden-to-hidden weights should have the same dimensions and values.
    assert(isequal(size(A), size(B)), '%s matrix size differs between %s %s and %s %s.', ...
        partName, labelA, mat2str(size(A)), labelB, mat2str(size(B)));
    maxDiff = maxAbsDiffLeading(A, B, size(A, 1), size(A, 2));
    assert(maxDiff <= tol, '%s weights differ between %s and %s. Max abs diff %.17g.', partName, labelA, labelB, maxDiff);
end

function maxDiff = maxAbsDiffLeading(A, B, nRows, nCols)
% Blocked comparison avoids allocating a second full hidden-to-hidden matrix.
    maxDiff = 0;
    blockRows = 512;
    for first = 1:blockRows:nRows
        last = min(first + blockRows - 1, nRows);
        d = max(abs(A(first:last, 1:nCols) - B(first:last, 1:nCols)), [], 'all');
        if d > maxDiff, maxDiff = d; end
    end
end
`;

const singleScript = String.raw`%% Check single-run network weight consistency
% Load top-level saved networks and verify that frozen initializer weights are
% identical where they should be shared across tasks.

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);
networkDir = fullfile(repoRoot, 'trained_networks');
files = dir(fullfile(networkDir, '*.mat'));
files = files(~startsWith({files.name}, '.'));
assert(~isempty(files), 'No single-run network artifacts found in trained_networks.');

records = loadWeightRecords(files);
compareAllRecords(records, 'Single-run network weight consistency');

${helpers}
`;

const seededScript = String.raw`%% Check seeded network weight consistency
% For each seed, verify that every saved seeded network uses the same frozen
% initializer samples. Hidden-to-hidden weights must match exactly. Encoder and
% decoder weights are compared over the greatest shared leading submatrix between
% networks, so tasks with different input or output counts can still be checked.
%
% The checker is written to avoid a full all-pairs comparison when possible. It
% first scans matrix sizes, then compares each encoder/decoder against the
% largest matrix that contains every other matrix as a leading submatrix. This is
% equivalent to all-pairs shared-submatrix checking for the deterministic
% initializer layout used here, but is much faster and uses less memory.

clear; clc;

repoRoot = locateRepoRoot();
addReleasePaths(repoRoot);

% These are the active seeded result folders. Add/remove entries here if you
% want to compare a different set of saved network collections.
SEED_COLLECTIONS = ["seeded", "seeded_state_random_derivative_ode45_lr_0p01"];

% Continuation training only changes learned biases, so the frozen weights should
% still match the same seed. Leave this true to include any continuation folders
% that are present locally.
AUTO_INCLUDE_DERIVATIVE_CONTINUATIONS = true;
OPTIONAL_SEED_COLLECTIONS = strings(0, 1);

TOLERANCE = 0;

if AUTO_INCLUDE_DERIVATIVE_CONTINUATIONS
    continuationDirs = dir(fullfile(repoRoot, 'trained_networks', 'seeded_state_random_derivative_ode45_continued_lr_0p001_extra_*'));
    continuationDirs = continuationDirs([continuationDirs.isdir]);
    SEED_COLLECTIONS = [SEED_COLLECTIONS, string({continuationDirs.name})];
end
SEED_COLLECTIONS = unique([SEED_COLLECTIONS(:); OPTIONAL_SEED_COLLECTIONS(:)], 'stable');
SEED_COLLECTIONS = SEED_COLLECTIONS(SEED_COLLECTIONS ~= "");

files = collectSeededNetworkFiles(repoRoot, SEED_COLLECTIONS);
assert(~isempty(files), 'No seeded network artifacts found in the selected trained_networks folders.');
seedValues = parseSeedValues(files);
uniqueSeeds = unique(seedValues, 'stable');

for s = 1:numel(uniqueSeeds)
    seed = uniqueSeeds(s);
    seedFiles = files(seedValues == seed);
    fprintf('\nChecking seed %s across %d networks.\n', seed, numel(seedFiles));
    compareSeedFileGroup(seedFiles, sprintf('Seeded network weight consistency for seed %s', seed), TOLERANCE);
end

fprintf('\nSeeded weight consistency passed for %d seeds across collections: %s.\n', numel(uniqueSeeds), strjoin(SEED_COLLECTIONS, ', '));

function files = collectSeededNetworkFiles(repoRoot, seedCollections)
    files = struct('folder', {}, 'name', {}, 'collection', {}, 'label', {});
    for c = 1:numel(seedCollections)
        collection = char(seedCollections(c));
        collectionRoot = fullfile(repoRoot, 'trained_networks', collection);
        if ~isfolder(collectionRoot)
            warning('Skipping missing seeded collection: %s', collectionRoot);
            continue
        end
        found = dir(fullfile(collectionRoot, '*', '*_seed_*_network.mat'));
        for i = 1:numel(found)
            [~, taskFolder] = fileparts(found(i).folder);
            files(end+1).folder = found(i).folder; %#ok<AGROW>
            files(end).name = found(i).name;
            files(end).collection = collection;
            files(end).label = sprintf('%s/%s/%s', collection, taskFolder, erase(found(i).name, '.mat'));
        end
    end
    if ~isempty(files)
        [~, order] = sort(string({files.label}));
        files = files(order);
    end
end

function seedValues = parseSeedValues(files)
    seedValues = strings(numel(files), 1);
    for i = 1:numel(files)
        token = regexp(files(i).name, '_seed_(\d+)_network\.mat$', 'tokens', 'once');
        assert(~isempty(token), 'Could not parse seed from %s.', files(i).name);
        seedValues(i) = string(token{1});
    end
end

function compareSeedFileGroup(files, titleText, tol)
    assert(numel(files) >= 2, '%s needs at least two networks to compare.', titleText);
    fprintf('%s\n', titleText);
    descriptors = scanWeightDescriptors(files);
    printDescriptorSummary(descriptors);
    compareExactPartByReference(files, descriptors, 'hidden', 'hidden-to-hidden', tol);
    compareSharedPartEfficient(files, descriptors, 'encoder', tol);
    compareSharedPartEfficient(files, descriptors, 'decoder', tol);
    fprintf('%s passed for %d networks.\n', titleText, numel(files));
end

function descriptors = scanWeightDescriptors(files)
    descriptors = struct('label', {}, 'encoderSize', {}, 'hiddenSize', {}, 'decoderSize', {});
    for i = 1:numel(files)
        modelPath = fullfile(files(i).folder, files(i).name);
        [net, label] = loadSavedNetwork(modelPath, files(i).label);
        sizes = extractWeightTripletSizes(net, label);
        descriptors(i).label = label; %#ok<AGROW>
        descriptors(i).encoderSize = sizes.encoder;
        descriptors(i).hiddenSize = sizes.hidden;
        descriptors(i).decoderSize = sizes.decoder;
        clear net
    end
end

function printDescriptorSummary(descriptors)
    for i = 1:numel(descriptors)
        fprintf('Found %-68s encoder %s | hidden %s | decoder %s\n', ...
            descriptors(i).label, mat2str(descriptors(i).encoderSize), ...
            mat2str(descriptors(i).hiddenSize), mat2str(descriptors(i).decoderSize));
    end
end

function compareExactPartByReference(files, descriptors, partName, displayName, tol)
    refIdx = 1;
    refLabel = descriptors(refIdx).label;
    refMatrix = loadWeightPartForFile(files(refIdx), partName);
    fprintf('Reference for %s: %s %s\n', displayName, refLabel, mat2str(size(refMatrix)));
    for i = 1:numel(files)
        if i == refIdx, continue; end
        candidate = loadWeightPartForFile(files(i), partName);
        maxDiff = compareExactMatrix(refMatrix, candidate, displayName, refLabel, descriptors(i).label, tol);
        fprintf('%s vs %s | %s %.3g\n', refLabel, descriptors(i).label, displayName, maxDiff);
        clear candidate
    end
    clear refMatrix
end

function compareSharedPartEfficient(files, descriptors, partName, tol)
    sizes = partSizes(descriptors, partName);
    [refIdx, hasDominant] = dominantReferenceIndex(sizes);
    if hasDominant
        refLabel = descriptors(refIdx).label;
        refMatrix = loadWeightPartForFile(files(refIdx), partName);
        fprintf('Reference for %s shared-submatrix checks: %s %s\n', partName, refLabel, mat2str(size(refMatrix)));
        for i = 1:numel(files)
            if i == refIdx, continue; end
            candidate = loadWeightPartForFile(files(i), partName);
            maxDiff = compareCommonMatrix(refMatrix, candidate, partName, refLabel, descriptors(i).label, tol);
            commonSize = [min(size(refMatrix, 1), size(candidate, 1)), min(size(refMatrix, 2), size(candidate, 2))];
            fprintf('%s vs %s | %s %.3g over %s\n', refLabel, descriptors(i).label, partName, maxDiff, mat2str(commonSize));
            clear candidate
        end
        clear refMatrix
    else
        warning('No single %s matrix contains all others as leading submatrices. Falling back to all-pairs shared-submatrix checks.', partName);
        compareSharedPartAllPairs(files, descriptors, partName, tol);
    end
end

function sizes = partSizes(descriptors, partName)
    sizes = zeros(numel(descriptors), 2);
    for i = 1:numel(descriptors)
        switch partName
            case 'encoder'
                sizes(i, :) = descriptors(i).encoderSize;
            case 'hidden'
                sizes(i, :) = descriptors(i).hiddenSize;
            case 'decoder'
                sizes(i, :) = descriptors(i).decoderSize;
            otherwise
                error('Unknown part name: %s', partName);
        end
    end
end

function [idx, hasDominant] = dominantReferenceIndex(sizes)
    maxRows = max(sizes(:, 1));
    maxCols = max(sizes(:, 2));
    candidates = find(sizes(:, 1) >= maxRows & sizes(:, 2) >= maxCols);
    hasDominant = ~isempty(candidates);
    if hasDominant
        [~, k] = min(prod(sizes(candidates, :), 2));
        idx = candidates(k);
    else
        [~, idx] = max(prod(sizes, 2));
    end
end

function compareSharedPartAllPairs(files, descriptors, partName, tol)
    for i = 1:numel(files)-1
        A = loadWeightPartForFile(files(i), partName);
        for j = i+1:numel(files)
            B = loadWeightPartForFile(files(j), partName);
            maxDiff = compareCommonMatrix(A, B, partName, descriptors(i).label, descriptors(j).label, tol);
            commonSize = [min(size(A, 1), size(B, 1)), min(size(A, 2), size(B, 2))];
            fprintf('%s vs %s | %s %.3g over %s\n', descriptors(i).label, descriptors(j).label, partName, maxDiff, mat2str(commonSize));
            clear B
        end
        clear A
    end
end

function matrix = loadWeightPartForFile(fileInfo, partName)
    modelPath = fullfile(fileInfo.folder, fileInfo.name);
    [net, label] = loadSavedNetwork(modelPath, fileInfo.label);
    matrix = extractWeightPart(net, label, partName);
    clear net
end

function [net, label] = loadSavedNetwork(modelPath, label)
    assert(isfile(modelPath), 'Missing saved network: %s', modelPath);
    S = load(modelPath);
    for candidate = ["net", "dlnet", "trainedNet", "netObj", "network", "bestNet"]
        name = char(candidate);
        if isfield(S, name) && isNetworkLike(S.(name))
            net = S.(name);
            return
        end
    end
    names = fieldnames(S);
    for k = 1:numel(names)
        if isNetworkLike(S.(names{k}))
            net = S.(names{k});
            return
        end
    end
    error('%s has no loadable network object.', modelPath);
end

function tf = isNetworkLike(v)
    c = lower(class(v));
    tf = contains(c, 'network') || strcmp(c, 'dlnetwork') || strcmp(c, 'seriesnetwork') || strcmp(c, 'dagnetwork');
end

function sizes = extractWeightTripletSizes(net, label)
    fcSizes = extractFullyConnectedWeightSizes(net);
    assert(size(fcSizes, 1) >= 3, '%s has %d fully connected weight matrices; expected at least 3.', label, size(fcSizes, 1));
    sizes.encoder = fcSizes(1, :);
    sizes.hidden = fcSizes(2, :);
    sizes.decoder = fcSizes(end, :);
end

function matrix = extractWeightPart(net, label, partName)
    [weightRefs, sourceType] = fullyConnectedWeightRefs(net);
    assert(numel(weightRefs) >= 3, '%s has %d fully connected weight matrices; expected at least 3.', label, numel(weightRefs));
    switch partName
        case 'encoder'
            idx = 1;
        case 'hidden'
            idx = 2;
        case 'decoder'
            idx = numel(weightRefs);
        otherwise
            error('Unknown part name: %s', partName);
    end
    switch sourceType
        case 'learnables'
            matrix = numericMatrix(net.Learnables.Value{weightRefs(idx)});
        case 'layers'
            matrix = numericMatrix(net.Layers(weightRefs(idx)).Weights);
        otherwise
            error('Unknown weight source type: %s', sourceType);
    end
end

function sizes = extractFullyConnectedWeightSizes(net)
    [weightRefs, sourceType] = fullyConnectedWeightRefs(net);
    sizes = zeros(numel(weightRefs), 2);
    for i = 1:numel(weightRefs)
        switch sourceType
            case 'learnables'
                W = net.Learnables.Value{weightRefs(i)};
            case 'layers'
                W = net.Layers(weightRefs(i)).Weights;
            otherwise
                error('Unknown weight source type: %s', sourceType);
        end
        sizes(i, :) = matrixSize(W);
    end
end

function [weightRefs, sourceType] = fullyConnectedWeightRefs(net)
    if isprop(net, 'Learnables')
        T = net.Learnables;
        params = string(T.Parameter);
        weightRefs = find(params == "Weights");
        if ~isempty(weightRefs)
            sourceType = 'learnables';
            return
        end
    end
    assert(isprop(net, 'Layers'), 'Network object has neither Learnables nor Layers.');
    layers = net.Layers;
    weightRefs = [];
    for i = 1:numel(layers)
        if isprop(layers(i), 'Weights') && ~isempty(layers(i).Weights)
            weightRefs(end+1) = i; %#ok<AGROW>
        end
    end
    sourceType = 'layers';
end

function sz = matrixSize(A)
    if isa(A, 'dlarray'), A = extractdata(A); end
    szRaw = size(A);
    assert(numel(szRaw) == 2 || all(szRaw(3:end) == 1), 'Expected a matrix weight tensor, got size %s.', mat2str(szRaw));
    sz = [szRaw(1), szRaw(2)];
end

function A = numericMatrix(A)
    if isa(A, 'dlarray'), A = extractdata(A); end
    try, A = gather(A); catch, end
    assert(isnumeric(A), 'Expected numeric weights, got %s.', class(A));
    assert(ismatrix(A), 'Expected a matrix weight tensor, got size %s.', mat2str(size(A)));
end

function maxDiff = compareCommonMatrix(A, B, partName, labelA, labelB, tol)
    nRows = min(size(A, 1), size(B, 1));
    nCols = min(size(A, 2), size(B, 2));
    assert(nRows > 0 && nCols > 0, '%s comparison between %s and %s has an empty common submatrix.', partName, labelA, labelB);
    maxDiff = maxAbsDiffLeading(A, B, nRows, nCols);
    assert(maxDiff <= tol, '%s weights differ between %s and %s over the greatest shared leading %d-by-%d submatrix. Max abs diff %.17g.', ...
        partName, labelA, labelB, nRows, nCols, maxDiff);
end

function maxDiff = compareExactMatrix(A, B, partName, labelA, labelB, tol)
    assert(isequal(size(A), size(B)), '%s matrix size differs between %s %s and %s %s.', ...
        partName, labelA, mat2str(size(A)), labelB, mat2str(size(B)));
    maxDiff = maxAbsDiffLeading(A, B, size(A, 1), size(A, 2));
    assert(maxDiff <= tol, '%s weights differ between %s and %s. Max abs diff %.17g.', partName, labelA, labelB, maxDiff);
end

function maxDiff = maxAbsDiffLeading(A, B, nRows, nCols)
    maxDiff = 0;
    blockRows = 512;
    for first = 1:blockRows:nRows
        last = min(first + blockRows - 1, nRows);
        blockDiff = abs(A(first:last, 1:nCols) - B(first:last, 1:nCols));
        d = max(blockDiff(:));
        if d > maxDiff, maxDiff = d; end
    end
    maxDiff = double(maxDiff);
end

function repoRoot = locateRepoRoot()
    starts = string.empty;
    scriptPath = mfilename('fullpath');
    if ~isempty(scriptPath)
        starts(end+1) = string(fileparts(scriptPath)); %#ok<AGROW>
    end
    stack = dbstack('-completenames');
    for k = 1:numel(stack)
        if isfield(stack(k), 'file') && ~isempty(stack(k).file)
            starts(end+1) = string(fileparts(stack(k).file)); %#ok<AGROW>
        end
    end
    starts(end+1) = string(pwd);
    for s = starts
        candidate = char(s);
        while ~isempty(candidate)
            if isfolder(fullfile(candidate, 'trained_networks')) && isfolder(fullfile(candidate, 'examples'))
                repoRoot = candidate;
                return
            end
            parent = fileparts(candidate);
            if strcmp(parent, candidate), break; end
            candidate = parent;
        end
    end
    error('Could not locate release root. Run from inside the release folder.');
end

function addReleasePaths(repoRoot)
    addpath(genpath(fullfile(repoRoot, 'src')));
    addpath(genpath(fullfile(repoRoot, 'examples')));
    addpath(genpath(fullfile(repoRoot, 'src', 'dynamical_systems')));
end
`;

makeMlx(path.join(repo, "tests", "check_single_network_weight_consistency.mlx"), singleScript);
makeMlx(path.join(repo, "tests", "check_seeded_network_weight_consistency.mlx"), seededScript);
console.log("Created weight-consistency live scripts.");
