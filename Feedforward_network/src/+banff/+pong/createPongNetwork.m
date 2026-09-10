function layers = createPongNetwork(width)
    layers = [
        featureInputLayer(5, "Name", "input")
        fullyConnectedLayer(width, "Name", "fc1", 'WeightsInitializer', @customWeights1_git, 'WeightLearnRateFactor', 0)
        tanhLayer("Name", "tanh1")
        fullyConnectedLayer(width, "Name", "fc2", 'WeightsInitializer', @customWeights2_git, 'WeightLearnRateFactor', 0)
        tanhLayer("Name", "tanh2")
        fullyConnectedLayer(3, "Name", "fc3", 'WeightsInitializer', @customWeights3_git, 'WeightLearnRateFactor', 0)
        softmaxLayer("Name", "softmax")
        classificationLayer];
end
