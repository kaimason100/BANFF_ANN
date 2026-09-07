function cfg = settings(task, overrides)
%SETTINGS One source of defaults for local and Slurm training.
if nargin < 2, overrides = struct(); end
cfg = struct('Width',16000,'MaxEpochs',1000,'LearnRate',0.01, ...
    'MiniBatchSize',500,'SplitSeed',42,'Plots',true, ...
    'ExecutionEnvironment','auto','NetworkSet','seeded_grouped');
switch char(task)
    case 'abalone', cfg.SplitSeed = 1;
    case 'Pong'
        cfg.MaxEpochs = 20000; cfg.LearnRate = 0.001; cfg.MiniBatchSize = 1000;
    case 'LQR_two_link_arm'
        cfg.MaxEpochs = 10000; cfg.MiniBatchSize = 5000;
end
names = fieldnames(overrides);
for k = 1:numel(names)
    assert(isfield(cfg,names{k}), 'Unknown training option: %s',names{k});
    cfg.(names{k}) = overrides.(names{k});
end
validateattributes(cfg.MaxEpochs,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.Width,{'numeric'},{'scalar','integer','positive','<=',16000});
validateattributes(cfg.LearnRate,{'numeric'},{'scalar','positive','finite'});
end
