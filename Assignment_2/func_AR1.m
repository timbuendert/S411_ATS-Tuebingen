function [AR1_r] = func_AR1(c, phi, yo, T)
% function to simulate an AR 1 process based on given parameters

epsilon = randn(T,1); % generate draws from SND as innovations
AR1_r = zeros(T,1); % initialize column vector for the realizations
for i = 1:T % generate realization for length T
    if i == 1 % need to include yo for t=1
        AR1_r(i,1) = c + phi * yo + epsilon(i);
    else
        AR1_r(i,1) = c + phi * AR1_r(i-1,1) + epsilon(i);
    end
end 
end