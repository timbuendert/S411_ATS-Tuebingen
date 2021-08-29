function [CLL_value] = lik_func_min(param, y)
% lik_func_min(param, y): same as @lik_func, but returns the negative value 
% of the total conditional log-likelihood function because these
% algorithms are minimization algorithms

% compute log likelihood contributions vector based on parameters and data
l_vector = lik_contrib(param, y);

% sum up l to compute value of the conditional log likelihood function and
% take the negative value of it
CLL_value = -sum(l_vector);
end