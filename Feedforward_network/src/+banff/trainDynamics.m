function trainDynamics(task,seeds,options)
%TRAINDYNAMICS One derivative-field trainer for local and Slurm entry points.
if nargin < 2, seeds = 0:9; end
if nargin < 3, options = struct(); end
names = {'Lorenz','MO0','MO5','MO7','MO13','Rikitake','SprottB','SprottC','SprottS'};
assert(any(strcmp(char(task),names)), 'Unknown dynamical system.');
dyn = readtable(fullfile(banff.root(),'examples','dynamical_systems','dynamics_list.xlsx'));
row = find(strcmp(string(dyn{:,1}),string(task)));
assert(isscalar(row),'The dynamics spreadsheet must contain exactly one matching row.');
validateattributes(seeds,{'numeric'},{'vector','integer','>=',0,'<=',9,'nonempty'});
for seed = seeds(:).'
    run_dynamical_system_state_random_derivative_ode45(row,seed,options);
end
end
