function drawCurvedArrow(startPos, endPos, curvature, color, lineWidth,arrowHeadSize)
    % DRAWCURVEDARROW Draws a curved arrow in MATLAB.
    % 
    % Inputs:
    %   startPos  - [x, y] start position
    %   endPos    - [x, y] end position
    %   curvature - Curvature factor (positive for counterclockwise, negative for clockwise)
    %   color     - Color of the arrow (e.g., 'k' for black)
    %   lineWidth - Width of the arrow line
    %
    % Example:
    %   banff.shared.drawCurvedArrow([0, 0], [1, 1], 0.5, 'r', 2);


    % Compute midpoint
    midPoint = (startPos + endPos) / 2;

    % Perpendicular direction for control point
    direction = endPos - startPos;
    normal = [-direction(2), direction(1)]; % Rotate 90 degrees
    normal = normal / norm(normal); % Normalize

    % Control point for Bézier curve
    controlPoint = midPoint + curvature * normal;

    % Generate Bézier curve points
    t = linspace(0, 1, 100);
    curveX = (1-t).^2 * startPos(1) + 2*(1-t).*t * controlPoint(1) + t.^2 * endPos(1);
    curveY = (1-t).^2 * startPos(2) + 2*(1-t).*t * controlPoint(2) + t.^2 * endPos(2);

    % Plot the curved line
    plot(curveX, curveY, 'Color', color, 'LineWidth', lineWidth);

    % Compute arrowhead direction
    arrowVec = [curveX(end) - curveX(end-1), curveY(end) - curveY(end-1)];
    arrowVec = arrowVec / norm(arrowVec); % Normalize

    % Define arrowhead size
    arrowSize = arrowHeadSize; % Scale with length
    leftTip = endPos - arrowSize * (arrowVec + [arrowVec(2), -arrowVec(1)] / 2);
    rightTip = endPos -arrowSize * (arrowVec - [arrowVec(2), -arrowVec(1)] / 2);

    % Draw arrowhead
    fill([endPos(1), leftTip(1), rightTip(1)], [endPos(2), leftTip(2), rightTip(2)], ...
        color, 'EdgeColor', color);

end
