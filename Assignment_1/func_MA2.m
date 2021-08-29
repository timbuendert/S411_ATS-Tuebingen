function [MA2_r] = func_MA2(mu, theta1, theta2, T)
% Function to generate realization of MA(2) process with given parameters

% generate draws from SND as innovations
epsilon = randn(T+2,1); % need T+2 because we need t = -2 and t = -1 for t = 0

% generate the realization by matrix computation
% integrate the lagged epsilons by sliding over the epsilon vector
% such that for t = 1, theta1 * epsilon(2 -> 0)+ theta2 * epsilon(1 -> -1) + epsilon(3 -> 1)
% such that for t = T, theta1 * epsilon(T+1 -> T-1)+ theta2 * epsilon(T -> T-2) + epsilon(T+2 -> T)

MA2_r = mu * ones(T,1) + theta1 * epsilon(2:T+1)+ theta2 * epsilon(1:T) + epsilon(3:T+2);

% alternative approach using a loop instead of matrix computations
%for i = 1:size(MA2_r,1)
%    MA2_r(i) = mu + phi1 * epsilon(i) + phi2 * epsilon(i+1) + epsilon(i+2);
%end