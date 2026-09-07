function replayPong(plotData,k,ax)
% replayPong  Animate a pong game from recorded plotData
%
%   banff.shared.replayPong(plotData) plays back the entire game stored in plotData,
%   showing the left "Opponent" paddle (blue) mirroring the ball and the
%   right "Network" paddle (green).  Ball is red.  A dashed mid-line and
%   labels indicate each side.

% Unpack
dt    = plotData.dt;
if isfield(plotData, 'paddleDisplayHeight')
    ph = plotData.paddleDisplayHeight;
else
    ph = 2 * plotData.paddleHeight / 1.5;
end
pw    = plotData.paddleWidth;
br    = plotData.ballRadius;
ballX = plotData.ballPosSeq(:,1);
ballY = plotData.ballPosSeq(:,2);
oppY  = plotData.oppPaddleYSeq;
netY  = plotData.netPaddleYSeq;
N     = numel(ballX);

% Paddle X positions
oppX = 0.05;
netX = 1 - oppX - pw;

% Playback loop
cla(ax);                % clear previous frame
hold(ax,'on');

% mid-court line and labels
plot(ax,[0.5 0.5],[0 1],'k--','LineWidth',1);
text(ax,0.25,0.95,'Opponent','Color','b', ...
    'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
text(ax,0.75,0.95,'Network','Color','g', ...
    'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');

% draw paddles
rectangle(ax,'Position',[oppX,   oppY(k)-ph/2, pw, ph], ...
    'FaceColor','b','EdgeColor','none');
rectangle(ax,'Position',[netX,   netY(k)-ph/2, pw, ph], ...
    'FaceColor','g','EdgeColor','none');

% draw ball
rectangle(ax,'Position',[ballX(k)-br, ballY(k)-br, 2*br, 2*br], ...
    'Curvature',[1,1],'FaceColor','r','EdgeColor','none');

hold(ax,'off');
%pause(dt);              % match original simulation speed
end
