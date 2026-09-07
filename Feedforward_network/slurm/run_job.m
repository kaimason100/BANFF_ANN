function run_job(task,seed,mode)
%RUN_JOB Thin headless adapter. No training or data-processing code lives here.
if nargin < 3, mode = 'train'; end
root = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root,'setup.m'));
set(groot,'DefaultFigureVisible','off');
task = char(task);
if strcmp(mode,'publication')
    assert(isempty(getenv('BANFF_LR')) && isempty(getenv('BANFF_EPOCHS')), ...
        'Publication mode uses the fixed schedule; remove BANFF_LR/BANFF_EPOCHS overrides.');
    if isempty(task)
        schedule = banff.publicationSchedule();
        tasks = unique(schedule.Task,'stable');
        validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<',10*numel(tasks)});
        task = char(tasks(floor(seed/10)+1));
        seed = mod(seed,10);
    end
    banff.runPublication(string(task),seed,false);
    return
end
ds = any(strcmp(task,{'Lorenz','MO0','MO5','MO7','MO13','Rikitake','SprottB','SprottC','SprottS'}));
options = struct();
if ~isempty(getenv('BANFF_LR')), options.LearnRate = str2double(getenv('BANFF_LR')); end
if ds
    if ~isempty(getenv('BANFF_EPOCHS')), options.TotalEpochs = str2double(getenv('BANFF_EPOCHS')); end
    if strcmp(mode,'continue')
        epochs = 50000; lr = 0.001;
        if isfield(options,'TotalEpochs'), epochs = options.TotalEpochs; end
        if isfield(options,'LearnRate'), lr = options.LearnRate; end
        options = banff.continuationOptions(epochs,lr);
    else
        assert(strcmp(mode,'train'),'Mode must be train or continue.');
    end
    banff.trainDynamics(task,seed,options);
else
    assert(strcmp(mode,'train'),'Continuation is supported only for derivative-field tasks.');
    options.Plots = false;
    if ~isempty(getenv('BANFF_EPOCHS')), options.MaxEpochs = str2double(getenv('BANFF_EPOCHS')); end
    banff.train(task,seed,options);
end
end
