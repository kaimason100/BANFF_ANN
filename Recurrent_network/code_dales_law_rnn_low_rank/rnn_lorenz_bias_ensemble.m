clearvars
close all
clc
%% Train and switch among the 20 Lorenz bias solutions used in Figure 3.
% All solutions share one fixed N=6000 network. Only their initial bias
% vectors differ. The first five solutions are switched every 400 time units
% in the long-timescale simulation. The rapid simulation selects among all 20
% solutions at each Euler step.
net_config
load network_configuration.mat

N_BIAS_SETS = 20;
TRAINING_ITERATIONS = 1e4;
LEARNING_RATE = 0.1;
MOMENTUM = 0.9; %#ok<NASGU> % Implemented in bias_gradient_descent.
INITIAL_BIAS_SD = 3;
LONG_SIMULATION_TIME = 2000;
LONG_SWITCH_INTERVAL = 400;
RAPID_DT = 1e-2;
RAPID_SIMULATION_TIME = 2000;

lorenzRow = find(strcmp(dynamics, 'Lorenz'), 1, 'first');
assert(~isempty(lorenzRow), 'Lorenz is missing from dynamics_list.xlsx.');
kLorenz = dim(lorenzRow);
assert(kLorenz == 3, 'The Figure 3 implementation expects a 3-D Lorenz system.');

% Compute the publication scaling and one fixed train/validation point set.
dt = 1e-2;
sampleTime = 1000;
sampleTimes = (0:round(sampleTime/dt))*dt;
y0 = ic(lorenzRow, 1:kLorenz)';
zeroSignal = zeros(size(sampleTimes));
[~, yScale] = ode45(@(t, y) int_dyn(y, 'Lorenz', sampleTimes, zeroSignal, t, 'Simulate'), ...
    [0 sampleTime], y0);
scal_mat = diag(max(abs(yScale), [], 1));
mean_mat = mean(yScale, 1)';

xTrain = 2*rand(kLorenz, nx) - 1;
xValidation = 2*rand(kLorenz, 100) - 1;
xTrainRaw = scal_mat*xTrain + mean_mat;
xValidationRaw = scal_mat*xValidation + mean_mat;
gTrain = pinv(scal_mat)*int_dyn(xTrainRaw, 'Lorenz', sampleTimes, ...
    zeros(1, nx), 0, 'Train') + xTrain;
gValidation = pinv(scal_mat)*int_dyn(xValidationRaw, 'Lorenz', sampleTimes, ...
    zeros(1, 100), 0, 'Train') + xValidation;

biasSets = zeros(N, N_BIAS_SETS);
trainingLoss = cell(N_BIAS_SETS, 1);
for biasIndex = 1:N_BIAS_SETS
    fprintf('Training Lorenz bias set %d/%d\n', biasIndex, N_BIAS_SETS);
    initialBias = INITIAL_BIAS_SD*randn(N, 1);
    [biasSets(:, biasIndex), ~, trainingLoss{biasIndex}] = ...
        bias_gradient_descent(xTrain, xValidation, zeros(1, nx), ...
        zeros(1, 100), win, eta, initialBias, phi, gTrain, gValidation, ...
        TRAINING_ITERATIONS, LEARNING_RATE, 1:N, []);
end

x0Norm = pinv(scal_mat)*(y0 - mean_mat);
longTimes = 0:dt:LONG_SIMULATION_TIME;
[tLong, xLong] = ode45(@(t, x) switchedLowRankDerivative(t, x, eta, phi, ...
    biasSets(:, 1:5), LONG_SWITCH_INTERVAL), longTimes, x0Norm);

rapidTimes = (0:RAPID_DT:RAPID_SIMULATION_TIME)';
xRapid = zeros(numel(rapidTimes), kLorenz);
xRapid(1, :) = x0Norm';
rapidBiasIndex = zeros(numel(rapidTimes) - 1, 1);
rapidStream = RandStream('mt19937ar', 'Seed', 1);
for step = 1:(numel(rapidTimes) - 1)
    rapidBiasIndex(step) = randi(rapidStream, N_BIAS_SETS);
    derivative = lowRankDerivative(xRapid(step, :)', eta, phi, ...
        biasSets(:, rapidBiasIndex(step)));
    xRapid(step + 1, :) = xRapid(step, :) + RAPID_DT*derivative';
end

save('lorenz_bias_ensemble_results.mat', 'biasSets', 'trainingLoss', 'eta', ...
    'phi', 'mean_mat', 'scal_mat', 'tLong', 'xLong', 'rapidTimes', ...
    'xRapid', 'rapidBiasIndex', 'TRAINING_ITERATIONS', 'LEARNING_RATE', ...
    'MOMENTUM', 'INITIAL_BIAS_SD', '-v7.3');

figure;
plot3(xLong(:, 1), xLong(:, 2), xLong(:, 3), 'k');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); grid on;
title('Lorenz output with bias switching every 400 time units');

figure;
plot(rapidTimes, xRapid);
xlabel('Time'); ylabel('Decoded state'); grid on;
title('Lorenz output with one randomly selected bias set per Euler step');

function dx = switchedLowRankDerivative(t, x, eta, phi, biasSets, switchInterval)
biasIndex = min(floor(t/switchInterval) + 1, size(biasSets, 2));
dx = lowRankDerivative(x, eta, phi, biasSets(:, biasIndex));
end

function dx = lowRankDerivative(x, eta, phi, bias)
dx = -x(:) + phi'*tanh(eta*x(:) + bias);
end
