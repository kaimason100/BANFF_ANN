function path = datasetPath(task)
%DATASETPATH Resolve the prepared data, including MATLAB's bundled Iris data.
task = char(task);
if strcmp(task,'iris')
    path = which('fisheriris.mat');
else
    if strcmp(task,'car_quality'), task = 'car'; end
    path = fullfile(banff.root(),'data',[task '_dataset.mat']);
end
assert(isfile(path),'Dataset not found for task %s.',task);
end
