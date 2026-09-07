function train(task, seeds, overrides)
%TRAIN Run the same seeded training on a workstation or a Slurm worker.
% Only Plots differs for headless jobs. Numerical options come from settings.
if nargin < 2, seeds = 0:9; end
if nargin < 3, overrides = struct(); end
validateattributes(seeds,{'numeric'},{'vector','integer','>=',0,'<=',9,'nonempty'});
assert(numel(unique(seeds)) == numel(seeds), 'Seeds must not contain duplicates.');
seeds = seeds(:).';
task = char(task);
cfg = banff.settings(task,overrides);
oldVisibility = get(groot,'DefaultFigureVisible');
cleanup = onCleanup(@() set(groot,'DefaultFigureVisible',oldVisibility)); %#ok<NASGU>
if cfg.Plots, set(groot,'DefaultFigureVisible','on'); else, set(groot,'DefaultFigureVisible','off'); end
if any(strcmp(task,{'iris','breast_cancer','car_quality','mushroom'}))
    banff.trainTabular(task,seeds,cfg);
elseif any(strcmp(task,{'abalone','toyota'}))
    banff.trainTabular(task,seeds,cfg);
elseif any(strcmp(task,{'mnist','mnist_fashion','kmnist','afro_mnist_ethiopic','afro_mnist_nko','afro_mnist_osmanya','afro_mnist_vai'}))
    banff.trainImages(task,seeds,cfg);
elseif strcmp(task,'Pong')
    for seed = seeds(:).', banff.trainPong(seed,cfg); end
elseif strcmp(task,'LQR_two_link_arm')
    for seed = seeds(:).', banff.trainMotor(seed,cfg); end
else
    error('Unknown feedforward task: %s. Use banff.trainDynamics for a dynamical system.',task);
end
end
