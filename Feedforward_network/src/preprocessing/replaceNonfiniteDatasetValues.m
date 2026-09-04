function X = replaceNonfiniteDatasetValues(X, datasetName, stageName)
%REPLACENONFINITEDATASETVALUES Report and replace non-finite data values.
%   X = REPLACENONFINITEDATASETVALUES(X, DATASETNAME, STAGENAME) preserves
%   the released zero-replacement behaviour while making every replacement
%   visible. The two documented categorical encodings receive informational
%   messages; all other replacements generate a warning.

if nargin < 2 || strlength(string(datasetName)) == 0
    datasetName = "unspecified dataset";
end
if nargin < 3 || strlength(string(stageName)) == 0
    stageName = "unspecified processing stage";
end

datasetName = string(datasetName);
stageName = string(stageName);
assert(isscalar(datasetName) && isscalar(stageName), ...
    'Dataset and processing-stage names must be scalar text values.');
datasetKey = lower(strrep(datasetName, "-", "_"));
nonfiniteMask = ~isfinite(X);

if any(datasetKey == ["car", "car_quality", "car_evaluation"]) && ...
        stageName == "loaded encoded matrix"
    X = reportCarEncoding(X, nonfiniteMask);
    return
end

if datasetKey == "mushroom" && stageName == "loaded encoded matrix"
    X = reportMushroomEncoding(X, nonfiniteMask);
    return
end

nReplaced = nnz(nonfiniteMask);
if nReplaced > 0
    warning('BANFF:UnexpectedDatasetImputation', ...
        ['%s: replacing %d unexpected non-finite value(s) with zero at ', ...
         '%s. Inspect the source data and preprocessing before using these results.'], ...
        char(datasetName), nReplaced, char(stageName));
    X(nonfiniteMask) = 0;
end
end

function X = reportCarEncoding(X, nonfiniteMask)
expectedRows = 1728;
expectedDoorCount = 432;
expectedPersonsCount = 576;

doorMask = false(size(X));
personsMask = false(size(X));
if ismatrix(X) && size(X, 2) >= 4
    doorMask(:, 3) = isnan(X(:, 3));
    personsMask(:, 4) = isnan(X(:, 4));
end
knownMask = doorMask | personsMask;
unknownMask = nonfiniteMask & ~knownMask;

nDoors = nnz(doorMask);
nPersons = nnz(personsMask);
matchesReleasedData = size(X, 1) == expectedRows && nDoors == expectedDoorCount && ...
    nPersons == expectedPersonsCount && nnz(unknownMask) == 0 && ...
    nnz(nonfiniteMask) == expectedDoorCount + expectedPersonsCount;

if matchesReleasedData
    fprintf(['[data handling] Car Evaluation: replacing %d doors=5more and ', ...
        '%d persons=more NaN category placeholders with reserved code 0 ', ...
        'before splitting and standardisation. This is categorical encoding, ', ...
        'not missing-data imputation.\n'], nDoors, nPersons);
else
    warning('BANFF:UnexpectedDatasetImputation', ...
        ['Car Evaluation encoding differs from the released dataset. Expected ', ...
         '%d doors=5more NaNs and %d persons=more NaNs in %d rows; observed ', ...
         '%d, %d and %d rows. All %d non-finite value(s) will be replaced ', ...
         'with zero; inspect the dataset before using these results.'], ...
        expectedDoorCount, expectedPersonsCount, expectedRows, nDoors, ...
        nPersons, size(X, 1), nnz(nonfiniteMask));
end

X(nonfiniteMask) = 0;
end

function X = reportMushroomEncoding(X, nonfiniteMask)
expectedRows = 8124;
expectedUnknownCount = 2480;
unknownCodeCount = 0;
if ismatrix(X) && size(X, 2) >= 12
    unknownCodeCount = nnz(X(:, 12) == 5);
end

matchesReleasedData = size(X, 1) == expectedRows && size(X, 2) == 23 && ...
    unknownCodeCount == expectedUnknownCount && ~any(nonfiniteMask(:));
if matchesReleasedData
    fprintf(['[data handling] Mushroom: %d unknown stalk-root entries are ', ...
        'already encoded as categorical code 5. No runtime imputation is ', ...
        'performed.\n'], unknownCodeCount);
else
    warning('BANFF:UnexpectedDatasetImputation', ...
        ['Mushroom encoding differs from the released dataset. Expected %d ', ...
         'rows, 23 columns, %d stalk-root code-5 entries and no non-finite ', ...
         'values; observed %d rows, %d columns, %d code-5 entries and %d ', ...
         'non-finite value(s). Any non-finite values will be replaced with zero.'], ...
        expectedRows, expectedUnknownCount, size(X, 1), size(X, 2), ...
        unknownCodeCount, nnz(nonfiniteMask));
end

X(nonfiniteMask) = 0;
end
