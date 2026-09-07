function [features, labels] = generateTrainingData(numSamples, ballSpeed)
    features = zeros(numSamples, 5);
    threshold = 0.1; centerPos = 0.5;
    labels = categorical(zeros(numSamples, 1), [0 1 2], ["Up", "Down", "Stay"]);
    for i = 1:numSamples
        ballX = rand; ballY = rand; paddleY = rand;
        theta = 2*pi*rand;
        ballVelX = ballSpeed * cos(theta); ballVelY = ballSpeed * sin(theta);
        features(i, :) = [ballX, ballY, ballVelX, ballVelY, paddleY];
        if ballX < 0.7
            if paddleY < centerPos - threshold, action = "Down";
            elseif paddleY > centerPos + threshold, action = "Up";
            else, action = "Stay"; end
        else
            if ballY < paddleY - threshold/2, action = "Up";
            elseif ballY > paddleY + threshold/2, action = "Down";
            else, action = "Stay"; end
        end
        labels(i) = categorical(action, ["Up", "Down", "Stay"]);
    end
end
