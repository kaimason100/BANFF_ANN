function dX = int_dyn(X, label, time, sig, t, flag)
%INT_DYN   Compute state derivatives for a variety of dynamical systems.
%   dX = INT_DYN(X, LABEL, TIME, SIG, T, FLAG) returns the time
%   derivatives of state(s) X for the system specified by LABEL.
%
%   Inputs:
%     X     — nStates×nSamples matrix of current state vectors.
%     LABEL — string name of the dynamical system to simulate:
%             'Integrator', 'Rossler', 'Lorenz', 'Thomas', 'SprottA'…'SprottS',
%             'Vanderpol', 'Pitchfork', 'Hopf Normal Form',
%             'Chua1'…'Chua6', 'Rikitake', 'Nose Hoover', 'Halvorsen',
%             'MO0'…'MO15'.
%     TIME  — 1×nSamples vector of time points corresponding to SIG.
%     SIG   — 1×nSamples external input signal.
%     T     — scalar or 1×m vector of query times at which to evaluate SIG.
%     FLAG  — 'Train' to use SIG directly, 'Simulate' to interpolate SIG.
%
%   Output:
%     dX    — nStates×m matrix of state derivatives at times T.
%
%   Example:
%     % Compute Lorenz derivatives at time tCurr
%     dX = int_dyn(X, 'Lorenz', tVec, extSignal, tCurr, 'Simulate');
%
%   Notes:
%     • RNG seed is not changed here; input-driven in other code.
%     • Ensure LABEL is one of the supported systems.
%     • All 3D systems are normalized by their TX time scale.

    % Determine drive input u depending on mode and signal
    if strcmp(flag, 'Train')           % training mode: use full signal
        u = sig;
    end

    if max(abs(sig)) > 0
        if strcmp(flag, 'Simulate')    % simulation mode: nearest‐neighbor interp
            u = interp1(time, sig, t, 'nearest');
        end
    else
        u = 0;                         % no input if signal is zero
    end

    %% 1D System: simple integrator
    if strcmp(label, 'Integrator')
        dX(1, :) = u;
        return;
    end

    %% 3D Chaotic Attractors & Other 3D Systems
    % Pre-extract state components for readability
    x = X(1, :);
    y = X(2, :);
    z = X(3, :);

    switch label

        case 'Rossler'
            % Rössler attractor parameters
            a = 0.2; b = 0.2; c = 5.7; TX = 10;
            dX(1, :) = -y - z;
            dX(2, :) = x + a * y;
            dX(3, :) = b + z .* (x - c);
            dX = dX / TX;

        case 'Lorenz'
            % Lorenz system parameters
            rho = 28; sigma = 10; beta = 8/3; TX = 40;
            dX(1, :) = sigma * (y - x);
            dX(2, :) = x .* (rho - z) - y;
            dX(3, :) = x .* y - beta * z;
            dX = dX / TX;

        case 'Thomas'
            % Thomas cyclically symmetric attractor
            b = 0.208186; TX = 1;
            dX(1, :) = sin(y) - b * x;
            dX(2, :) = sin(z) - b * y;
            dX(3, :) = sin(x) - b * z;
            dX = dX / TX;

        case {'SprottA','SprottB','SprottC','SprottD','SprottE','SprottF', ...
              'SprottG','SprottH','SprottI','SprottJ','SprottK','SprottL', ...
              'SprottM','SprottN','SprottO','SprottP','SprottQ','SprottR','SprottS'}
            % Family of simple 3D "Sprott" systems
            TX = 10;
            switch label
                case 'SprottA'
                    dX(1,:) = y;
                    dX(2,:) = -x + y .* z;
                    dX(3,:) = 1 - y.^2;
                case 'SprottB'
                    dX(1,:) = y .* z;
                    dX(2,:) = x - y;
                    dX(3,:) = 1 - x .* y;
                case 'SprottC'
                    dX(1,:) = y .* z;
                    dX(2,:) = x - y;
                    dX(3,:) = 1 - x.^2;
                case 'SprottD'
                    dX(1,:) = -y;
                    dX(2,:) = x + z;
                    dX(3,:) = x .* z + 3 * y.^2;
                case 'SprottE'
                    dX(1,:) = y .* z;
                    dX(2,:) = x.^2 - y;
                    dX(3,:) = 1 - 4 * x;
                case 'SprottF'
                    dX(1,:) = y + z;
                    dX(2,:) = -x + 0.5 * y;
                    dX(3,:) = x.^2 - z;
                case 'SprottG'
                    dX(1,:) = 0.4 * x + z;
                    dX(2,:) = x .* z - y;
                    dX(3,:) = -x + y;
                case 'SprottH'
                    dX(1,:) = -y + z.^2;
                    dX(2,:) = x + 0.5 * y;
                    dX(3,:) = x - z;
                case 'SprottI'
                    % Note: two variants exist—ensure correct formula.
                    dX(1,:) = -0.2 * y;
                    dX(2,:) = x + z;
                    dX(3,:) = x + y.^2 - z;
                case 'SprottJ'
                    dX(1,:) = 2 * z;
                    dX(2,:) = -2 * y + z;
                    dX(3,:) = -x + y + y.^2;
                case 'SprottK'
                    dX(1,:) = x .* y - z;
                    dX(2,:) = x - y;
                    dX(3,:) = x + 0.3 * z;
                case 'SprottL'
                    TX = 1;
                    dX(1,:) = y + 3.9 * z;
                    dX(2,:) = 0.9 * x.^2 - y;
                    dX(3,:) = 1 - x;
                case 'SprottM'
                    dX(1,:) = -z;
                    dX(2,:) = -x.^2 - y;
                    dX(3,:) = 1.7 + 1.7 * x + y;
                case 'SprottN'
                    dX(1,:) = -2 * y;
                    dX(2,:) = x + z.^2;
                    dX(3,:) = 1 + y - 2 * z;
                case 'SprottO'
                    dX(1,:) = y;
                    dX(2,:) = x - z;
                    dX(3,:) = x + x .* z + 2.7 * y;
                case 'SprottP'
                    dX(1,:) = 2.7 * y + z;
                    dX(2,:) = -x + y.^2;
                    dX(3,:) = x + y;
                case 'SprottQ'
                    dX(1,:) = -z;
                    dX(2,:) = x - y;
                    dX(3,:) = 3.1 * x + y.^2 + 0.5 * z;
                case 'SprottR'
                    dX(1,:) = 0.9 - y;
                    dX(2,:) = 0.4 + z;
                    dX(3,:) = x .* y - z;
                case 'SprottS'
                    dX(1,:) = -x - 4 * y;
                    dX(2,:) = x + z.^2;
                    dX(3,:) = 1 + x;
            end
            dX = dX / TX;

        case 'Vanderpol'
            % Van der Pol oscillator
            mu = 5; TX = 10;
            dX(1,:) = mu * (x - (x.^3)/3 - y);
            dX(2,:) = x / mu;
            dX = dX / TX;

        case 'Pitchfork'
            % Pitchfork bifurcation normal form
            TX = 1;
            dX(1,:) = 0.5 * x - x.^3;
            dX = dX / TX;

        case 'Hopf Normal Form'
            % Hopf bifurcation in normal form
            beta = 0.5; sigma = -1; TX = 10;
            dX(1,:) = beta * x - y + sigma * x .* (x.^2 + y.^2);
            dX(2,:) = x + beta * y + sigma * y .* (x.^2 + y.^2);
            dX = dX / TX;

        case {'Chua1','Chua2','Chua3','Chua4','Chua5','Chua6'}
            % Variants of Chua's circuit
            TX = 10;
            switch label
                case 'Chua1'
                    dX(1,:) = 0.3 * y + x - x.^3;
                case 'Chua2'
                    dX(1,:) = 0.2 * y - x + 2 * tanh(x);
                case 'Chua3'
                    dX(1,:) = 0.2 * y + x - x .* abs(x);
                case 'Chua4'
                    dX(1,:) = 0.2 * y - x - 2 * sin(x);
                case 'Chua5'
                    dX(1,:) = 0.2 * y - 0.3 * x + sign(x);
                case 'Chua6'
                    dX(1,:) = 0.2 * y - x + 2 * atan(x);
            end
            dX(2,:) = x + z;
            dX(3,:) = -y;
            dX = dX / TX;

        case 'Rikitake'
            % Rikitake two-disk dynamo
            mu = 1; alpha = 1; TX = 10;
            dX(1,:) = -mu * x + y .* z;
            dX(2,:) = -mu * y + x .* (z - alpha);
            dX(3,:) = 1 - x .* y;
            dX = dX / TX;

        case 'Nose Hoover'
            % Nosé–Hoover thermostat
            TX = 10;
            dX(1,:) = y;
            dX(2,:) = y .* z - x;
            dX(3,:) = 1 - y.^2;
            dX = dX / TX;

        case 'Halvorsen'
            % Halvorsen attractor (note: parameter a must be defined)
            TX = 10;
            a = 1;  % <-- ensure a is set appropriately
            dX(1,:) = -a * x - 4 * y - 4 * z - y.^2;
            dX(2,:) = -a * y - 4 * z - 4 * x - z.^2;
            dX(3,:) = a * z - 4 * x - 4 * y - x.^2;
            dX = dX / TX;

        case {'MO0','MO1','MO2','MO3','MO4','MO5','MO6','MO7','MO8','MO9','MO10','MO11','MO12','MO13','MO14','MO15'}
            % Family of "MO" systems (various nonlinear oscillators)
            TX = 10; a = 0.6; b = 1;
            switch label
                case 'MO0'
                    g = abs(x) - 1;
                case 'MO1'
                    g = 1 - 6 * max(x, 0);
                case 'MO2'
                    g = sign(x) - x;
                case 'MO3'
                    a = 1; g = 1.1 * (x.^2 - 1);
                case 'MO4'
                    g = x .* (x - 1);
                case 'MO5'
                    g = x .* (1 - x.^2);
                case 'MO6'
                    g = (x.^2) .* (1 - x);
                case 'MO7'
                    g = (x.^2) .* (1 - x.^2);
                case 'MO8'
                    g = x .* (x.^4 - 1);
                case 'MO9'
                    g = (x.^3) .* (1 - x);
                case 'MO10'
                    g = (x.^2) .* (1 - x.^3);
                case 'MO11'
                    a = 1; g = 5 - exp(x);
                case 'MO12'
                    a = 1; g = 7 - 8 * tanh(x);
                case 'MO13'
                    g = 6 * tanh(x) - 3 * x;
                case 'MO14'
                    g = 6 * atan(x) - x;
                case 'MO15'
                    g = x - 0.5 * sinh(x);
            end
            dX(1,:) = y;
            dX(2,:) = z;
            dX(3,:) = g - a * z - b * y;
            dX = dX / TX;

        otherwise
            error('Unknown system label "%s".', label);
    end
end
