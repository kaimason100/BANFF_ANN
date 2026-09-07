function [A, B, K, dt, T, t, l1, l2] = lqrArmSetup()
    A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    B = [0 0; 0 0; 1 0; 0 1];
    Q = diag([100, 100, 1, 1]); R = diag([0.1, 0.1]);
    K = lqr(A, B, Q, R);
    dt = 0.01; T = 5; t = 0:dt:T; l1 = 1.0; l2 = 0.8;
end
