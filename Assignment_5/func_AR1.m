function [AR1_r] = func_AR1(c, phi, yo, trend, T)
% function to simulate an AR 1 process based on given parameters

epsilon = randn(T,1); % generate draws from SND as innovations
AR1_r = zeros(T+1,1); % initialize column vector for the realizations
AR1_r(1) = yo;
for i = 2:T+1 % generate realization for length of T
    AR1_r(i,1) = c + phi * AR1_r(i-1,1) + trend * (i-1) + epsilon(i-1) ;
end