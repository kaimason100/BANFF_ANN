function stop = progress(info, maxEpochs, tag)
%PROGRESS Display loss and elapsed/remaining time without consuming any RNG.
persistent started lastEpoch
stop = false;
if strcmp(info.State,'start')
    started = tic; lastEpoch = -1;
elseif strcmp(info.State,'iteration') && ...
        (info.Epoch == 1 || mod(info.Epoch,10) == 0 || info.Epoch == maxEpochs) && info.Epoch ~= lastEpoch
    lastEpoch = info.Epoch;
    elapsed = toc(started);
    eta = elapsed * (maxEpochs / max(1,info.Epoch) - 1);
    fprintf('[%s] epoch %d/%d | %.1f%% | loss %.6g | elapsed %.0f s | ETA %.0f s\n', ...
        tag,info.Epoch,maxEpochs,100*info.Epoch/maxEpochs,info.TrainingLoss,elapsed,eta);
elseif strcmp(info.State,'done')
    fprintf('[%s] finished | elapsed %.0f s\n',tag,toc(started));
end
end
