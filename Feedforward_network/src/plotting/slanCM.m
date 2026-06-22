function cmap = slanCM(name, n)
% slanCM  Local colormap compatibility for publication scripts.
%
% The original publication scripts used slanCM colormap names. This lightweight
% local implementation keeps the scripts self-contained when the third-party
% colormap package is not installed.

if nargin < 2 || isempty(n)
    n = 256;
end
requestedName = string(name);
name = lower(requestedName);

installedCmap = tryInstalledSlanCM(requestedName, n);
if ~isempty(installedCmap)
    cmap = installedCmap;
    fprintf('slanCM: requested "%s", using installed slanCM with %d colours.\n', requestedName, n);
    return;
end

resolvedName = name;
switch name
    case "iceburn"
        anchors = [
            0.0000 0.0000 0.0000
            0.1059 0.2353 0.4314
            0.2667 0.5725 0.7647
            0.9294 0.9294 0.9294
            0.8627 0.4314 0.2039
            0.5020 0.0784 0.0784
            0.0000 0.0000 0.0000
        ];
    case {"green", "greens"}
        anchors = [
            0.9686 0.9882 0.9608
            0.8980 0.9608 0.8784
            0.7804 0.9137 0.7529
            0.6314 0.8510 0.6078
            0.4549 0.7686 0.4627
            0.2549 0.6706 0.3647
            0.1373 0.5451 0.2706
            0.0000 0.4275 0.1725
            0.0000 0.2667 0.1059
        ];
    case "pink"
        cmap = pink(n);
        fprintf('slanCM: requested "%s", using "%s" with %d colours.\n', requestedName, resolvedName, n);
        return;
    otherwise
        if isMatlabColormap(name)
            cmap = feval(char(name), n);
            fprintf('slanCM: requested "%s", using MATLAB colormap "%s" with %d colours.\n', requestedName, name, n);
            return;
        end
        error('Unknown slanCM colormap "%s". Install the full slanCM package, add it to the MATLAB path ahead of this project, or add this map to src/plotting/slanCM.m.', requestedName);
end
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'linear');
cmap = min(max(cmap, 0), 1);
fprintf('slanCM: requested "%s", using "%s" with %d colours.\n', requestedName, resolvedName, n);
end

function cmap = tryInstalledSlanCM(name, n)
cmap = [];
currentDir = fileparts(mfilename('fullpath'));
pathParts = strsplit(path, pathsep);
removeMask = false(size(pathParts));
for i = 1:numel(pathParts)
    p = pathParts{i};
    if isempty(p), continue; end
    if strcmpi(p, currentDir) || endsWith(strrep(p, '\', '/'), '/src/plotting')
        removeMask(i) = true;
    end
end
if ~any(removeMask)
    return;
end
originalPath = path;
cleanupObj = onCleanup(@() path(originalPath)); %#ok<NASGU>
rmpath(pathParts{removeMask});
try
    if exist('slanCM', 'file') == 2
        cmap = slanCM(char(name), n);
    end
catch
    cmap = [];
end
end

function tf = isMatlabColormap(name)
tf = any(strcmp(name, ["parula", "turbo", "hsv", "hot", "cool", "spring", "summer", "autumn", "winter", "gray", "bone", "copper", "pink", "jet", "lines", "colorcube", "prism", "flag", "white"]));
end
