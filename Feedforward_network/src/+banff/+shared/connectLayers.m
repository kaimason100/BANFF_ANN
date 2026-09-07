function connectLayers(fromPositions, toPositions)
    for i = 1:size(fromPositions,1)
        for j = 1:size(toPositions,1)
            plot([fromPositions(i,1), toPositions(j,1)], ...
                 [fromPositions(i,2), toPositions(j,2)], ...
                 'Color','k', 'LineWidth',1);
        end
    end
end

%% Function to draw shaded box around layer
