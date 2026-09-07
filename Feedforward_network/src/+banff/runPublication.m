function runPublication(tasks,seeds,showPlots)
%RUNPUBLICATION Execute the fixed publication schedule, including continuation.
% Each task/seed resets exactly the same numerical RNGs whether called alone,
% in this local batch, or in a Slurm array. Continuation follows its source.
schedule = banff.publicationSchedule();
if nargin < 1 || isempty(tasks), tasks = unique(schedule.Task,'stable'); end
if nargin < 2, seeds = 0:9; end
if nargin < 3, showPlots = true; end
assert(all(ismember(string(tasks),schedule.Task)), 'Unknown publication task.');
validateattributes(seeds,{'numeric'},{'vector','integer','>=',0,'<=',9,'nonempty'});
selected = ismember(schedule.Task,string(tasks)) & ismember(schedule.Seed,seeds);
schedule = schedule(selected,:);
oldVisibility = get(groot,'DefaultFigureVisible');
cleanup = onCleanup(@() set(groot,'DefaultFigureVisible',oldVisibility)); %#ok<NASGU>
if showPlots, set(groot,'DefaultFigureVisible','on'); else, set(groot,'DefaultFigureVisible','off'); end
disp(schedule);
for k = 1:height(schedule)
    row = schedule(k,:);
    fprintf('[publication %d/%d] %s seed %d | %s | %d epochs | learning rate %.6g\n', ...
        k,height(schedule),row.Task,row.Seed,row.Stage,row.Epochs,row.LearnRate);
    options = struct('LearnRate',row.LearnRate);
    if row.DerivativeField
        options.TotalEpochs = row.Epochs;
        options.OutputNetworkSet = 'seeded_state_random_derivative_ode45_lr_0p01';
        if row.Stage == "continuation"
            options = banff.continuationOptions(row.Epochs,row.LearnRate);
        end
        banff.trainDynamics(row.Task,row.Seed,options);
    else
        options.MaxEpochs = row.Epochs; options.Plots = showPlots;
        banff.train(row.Task,row.Seed,options);
    end
end
end
