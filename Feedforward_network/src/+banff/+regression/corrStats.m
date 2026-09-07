function [r, p] = corrStats(yPred, yTrue)
    valid = isfinite(yPred) & isfinite(yTrue);
    yPred = yPred(valid); yTrue = yTrue(valid);
    if numel(yPred) < 3 || std(yPred) == 0 || std(yTrue) == 0
        r = NaN; p = NaN; return
    end
    [R, P] = corrcoef(yPred, yTrue); r = R(1, 2); p = P(1, 2);
end
