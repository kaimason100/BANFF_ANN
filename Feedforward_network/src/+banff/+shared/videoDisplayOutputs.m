function outputs = videoDisplayOutputs(plotData, taskIdx)
outputs = plotData.output;
if taskIdx == 2
    if isfield(plotData, 'actions')
        actions = string(plotData.actions(:));
        outputs = zeros(numel(actions), 3);
        outputs(:, 1) = actions == "Down";
        outputs(:, 2) = actions == "Up";
        outputs(:, 3) = actions == "Stay";
    elseif size(outputs, 2) == 1
        outputs = [outputs(:), zeros(numel(outputs), 2)];
    end
elseif taskIdx == 3 && size(outputs, 2) == 2
    outputs = [outputs, zeros(size(outputs, 1), 1)];
end
assert(size(outputs, 2) >= 3, 'Video display outputs must have at least three columns after task-specific formatting.');
outputs = double(outputs(:, 1:3));
end
