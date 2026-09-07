function value = plotOption(enabled)
%PLOTOPTION Only presentation differs between local and Slurm runs.
if enabled, value = 'training-progress'; else, value = 'none'; end
end
