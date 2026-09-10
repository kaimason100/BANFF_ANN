function assertDerivativeModelIdentity(S, task, seed, basePath, isContinuation)
%ASSERTDERIVATIVEMODELIDENTITY Reject mislabelled models and stale continuations.
% Legacy continuations lack source hashes; warn rather than imply verification.
if nargin < 4, basePath = ''; end
if nargin < 5, isContinuation = false; end
assert(isfield(S,'metadata'),'banff:ModelIdentity','Missing DS model metadata.');
M = S.metadata;
assert(isfield(M,'task') && isequal(string(M.task),string(task)), ...
    'banff:ModelIdentity','Requested %s seed %d, but saved task metadata differs.',task,seed);
assert(isfield(M,'seed') && isequal(double(M.seed),double(seed)), ...
    'banff:ModelIdentity','%s saved seed does not match requested seed %d.',task,seed);
assert(isfield(M,'trainingSamplingMode') && isequal(string(M.trainingSamplingMode),"state-random-derivative-field"), ...
    'banff:ModelIdentity','%s seed %d is not a derivative-field model.',task,seed);
if ~isContinuation, return; end
[taskFolder,~,~] = fileparts(char(basePath));
[setFolder,~,~] = fileparts(taskFolder);
[~,baseSet] = fileparts(setFolder);
assert(isfield(M,'sourceNetworkSet') && isequal(string(M.sourceNetworkSet),string(baseSet)), ...
    'banff:ModelIdentity','%s seed %d continuation belongs to a different base network set.',task,seed);
if isfield(M,'sourceMetadata') && isstruct(M.sourceMetadata) && ~isempty(fieldnames(M.sourceMetadata))
    banff.assertDerivativeModelIdentity(struct('metadata',M.sourceMetadata),task,seed);
end
if isfield(M,'sourceNetworkSHA256') && strlength(string(M.sourceNetworkSHA256)) > 0
    assert(isfile(basePath),'banff:ModelIdentity','Cannot verify continuation: base network is missing: %s',basePath);
    assert(strcmpi(banff.fileSHA256(basePath),char(M.sourceNetworkSHA256)), ...
        'banff:ModelIdentity','%s seed %d continuation source differs from the current base file.',task,seed);
else
    warning('banff:UnverifiedContinuationSource', ...
        '%s seed %d: legacy continuation has no source hash. Task, seed and base set match, but exact source-file provenance is unverified.',task,seed);
end
end
