function path = networkPath(repoRoot,networkSet,task,seed)
%NETWORKPATH Keep new grouped-split models separate from historical results.
folder = fullfile(repoRoot,'trained_networks',char(networkSet),char(task));
if ~isfolder(folder), mkdir(folder); end
path = fullfile(folder,sprintf('%s_seed_%03d_network.mat',char(task),seed));
end
