function nFrames = taskVideoFrameCount(plotData, taskIdx)
if taskIdx == 2
    nFrames = min([size(plotData.ballPosSeq, 1), numel(plotData.oppPaddleYSeq), numel(plotData.netPaddleYSeq)]);
elseif taskIdx == 3
    nFrames = min([size(plotData.netTrajectory, 2), size(plotData.lqrTrajectory, 2)]);
else
    nFrames = inf;
end
end
