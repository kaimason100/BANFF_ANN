function shadedBox(x, yCoords, color)
    yMin = min(yCoords) - 0.4;
    yMax = max(yCoords) + 0.4;
    rectangle('Position',[x-0.5, yMin, 1, yMax - yMin], 'FaceColor', color, 'EdgeColor', 'none','FaceAlpha',0.2);
end

%% Function to draw bias vector
