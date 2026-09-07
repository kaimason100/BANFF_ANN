function repoRoot = root()
%ROOT Resolve the package from this function, independent of current folder.
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
assert(isfolder(fullfile(repoRoot,'data')) && isfolder(fullfile(repoRoot,'examples')), ...
    'Incomplete Feedforward_network package.');
end
