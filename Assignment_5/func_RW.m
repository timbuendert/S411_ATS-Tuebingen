function [RW_r] = func_RW(y0, T)
% function to simulate Random Walk process (AR(1) with unit root) based 
% on given parameters

epsilon = randn(T+1,1); % generate draws from SND as innovations
RW_r = zeros(T+1,1); % initialize column vector for the realizations
RW_r(1,1) = y0;

for i = 2:T+1 % generate realization for length of T
    RW_r(i,1) = RW_r(i-1,1) + epsilon(i);
end 
end