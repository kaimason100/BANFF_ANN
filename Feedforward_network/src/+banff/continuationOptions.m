function options = continuationOptions(epochs,learnRate)
%CONTINUATIONOPTIONS Shared source, output naming and optimiser continuation.
if nargin < 1, epochs = 50000; end
if nargin < 2, learnRate = 0.001; end
validateattributes(epochs,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(learnRate,{'numeric'},{'scalar','positive','finite'});
label = strrep(sprintf('%.15g',learnRate),'.','p');
label = strrep(label,'-','m');
options = struct('TotalEpochs',epochs,'LearnRate',learnRate, ...
    'ContinueFromSource',true, ...
    'SourceNetworkSet','seeded_state_random_derivative_ode45_lr_0p01', ...
    'OutputNetworkSet',sprintf('seeded_state_random_derivative_ode45_continued_lr_%s_extra_%d',label,epochs), ...
    'ContinuationLabel',sprintf('continued_lr_%s_extra_%d',label,epochs));
end
