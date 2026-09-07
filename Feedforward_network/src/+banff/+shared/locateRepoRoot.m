function repoRoot = locateRepoRoot()
    starts = string.empty;
    scriptPath = mfilename('fullpath');
    if ~isempty(scriptPath)
        starts(end+1) = string(fileparts(scriptPath)); %#ok<AGROW>
    end
    stack = dbstack('-completenames');
    for k = 1:numel(stack)
        if isfield(stack(k), 'file') && ~isempty(stack(k).file)
            starts(end+1) = string(fileparts(stack(k).file)); %#ok<AGROW>
        end
    end
    starts(end+1) = string(pwd);
    for s = starts
        candidate = char(s);
        while ~isempty(candidate)
            if isfolder(fullfile(candidate, 'examples')) && isfolder(fullfile(candidate, 'src'))
                repoRoot = candidate;
                return
            end
            nestedCandidate = fullfile(candidate, 'Feedforward_network');
            if isfolder(fullfile(nestedCandidate, 'examples')) && isfolder(fullfile(nestedCandidate, 'src'))
                repoRoot = nestedCandidate;
                return
            end
            parent = fileparts(candidate);
            if strcmp(parent, candidate), break; end
            candidate = parent;
        end
    end
    error('Could not locate release root. Run from inside the release folder.');
end
