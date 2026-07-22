clearvars
close all
clc
tic
%% Train the Dale-constrained recurrent BANFF network for Figure S4.
net_config_dales_law
epsilon = 0.02;
no = 3*10^4;

load network_configuration_dales_law.mat
k0 = size(phi, 2);

% Figure S4 displays Rossler, Lorenz, Thomas, Rikitake, MO0 and SprottB.
% These indices refer to data rows after readtable removes the header row.
SYSTEM_ROWS = [5 6 7 33 35 9];
assert(all(SYSTEM_ROWS <= nsys), 'A requested system row exceeds the dynamics table.');
expectedSystems = ["Rossler"; "Lorenz"; "Thomas"; "Rikitake"; "MO0"; "SprottB"];
selectedSystems = strtrim(string(dynamics(SYSTEM_ROWS)));
assert(isequal(selectedSystems(:), expectedSystems), ...
    'SYSTEM_ROWS no longer selects the six systems displayed in Figure S4.');

for learn_sys = SYSTEM_ROWS
    x = 2*rand(k0, nx) - 1;
    xv = 2*rand(k0, 100) - 1;
    c = 2*rand(1, nx) - 1;
    cv = 2*rand(1, 100) - 1;

    dt = 1e-2;
    k = dim(learn_sys);
    bias0 = randn(N, 1);
    y0 = ic(learn_sys, 1:k)';
    sampleTime = 1000;
    sampleTimes = (0:round(sampleTime/dt))*dt;
    sampleSignal = cos(2*pi*sampleTimes'*40*rand(1, 20)/1000)*rand(20, 1);
    sampleSignal(abs(sampleSignal) > 1) = 0;

    label = dynamics{learn_sys};
    [~, yScale] = ode45(@(t, y) int_dyn(y, label, sampleTimes, ...
        input(learn_sys, 1)*sampleSignal, input(learn_sys, 1)*t, 'Simulate'), ...
        [0 sampleTime], y0);

    xk = x(1:k, :);
    xvk = xv(1:k, :);
    if rescale(learn_sys) == 1
        scal_mat = diag(max(abs(yScale), [], 1));
        mean_mat = mean(yScale, 1)';
    else
        scal_mat = eye(k);
        mean_mat = zeros(k, 1);
    end
    xscaled = scal_mat*xk + mean_mat;
    xscaledv = scal_mat*xvk + mean_mat;
    gx = pinv(scal_mat)*int_dyn(xscaled, label, sampleTimes, ...
        input(learn_sys, 1)*c, 0, 'Train') + xk;
    gxv = pinv(scal_mat)*int_dyn(xscaledv, label, sampleTimes, ...
        input(learn_sys, 1)*cv, 0, 'Train') + xvk;

    if k < k0
        i1 = find(sum(abs(eta(:, (k + 1):end)), 2) == 0);
        i2 = find(sum(abs(eta(:, (k + 1):end)), 2) > 0);
    else
        i1 = 1:N;
        i2 = [];
    end

    [bias, xhat, store, storeb] = bias_gradient_descent(xk, xvk, ...
        input(learn_sys, 1)*c, input(learn_sys, 1)*cv, win, eta, bias0, ...
        phi, gx, gxv, no, epsilon, i1, i2); %#ok<ASGLU>

    outputFile = sprintf('%s_dales_law_results.mat', matlab.lang.makeValidName(label));
    save(outputFile, 'bias', 'eta', 'phi', 'store', 'storeb', ...
        'mean_mat', 'scal_mat', 'epsilon', 'no', '-v7.3');
    close all
end
