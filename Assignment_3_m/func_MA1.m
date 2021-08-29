function [MA1_r] = func_MA1(mu, theta, sigma2, T)
% func_MA1(mu, theta, sigma2, T): function to simulate a realization (MA1_r)
% of MA(1) process with given parameters (μ, θ, σ^2) for 
% length T (t = 1, ..., T)

% general MA(1): y_t = μ + θ * ɛ_t-1 + ɛ_t

% generate draws from SND (as ɛ_t ~ N(0,σ^2)) for t = 0, ..., T -> T+1 draws
epsilon = normrnd(0, sqrt(sigma2), [T+1,1]); % take sqrt because normrnd(μ, s.d.)

% generate the realization by matrix computation: integrate the lagged 
% epsilons by sliding over the epsilon vector
% - such that for t = 1: y_1 = μ + θ * ɛ_0 [-> epsilon(1)] + ɛ_1 [-> epsilon(2)]
% ...
% - such that for t = T: y_T = μ + θ * ɛ_T-1 [-> epsilon(T)] + ɛ_T [-> epsilon(T+1)]

MA1_r = mu * ones(T,1) + theta * epsilon(1:T)+ epsilon(2:T+1);
end