function [CLL_value] = lik_func(param, y)
% lik_func(param, y): function to first, determine the vector of individual
% log-likelihood contributions based on parameters (param) and data
% (y), and, second, summing up all individual contributions to compute and 
% return the value of the total conditional log-likelihood function
% (CLL_value)

% compute log likelihood contributions vector based on parameters and data
l_vector = lik_contrib(param, y);

% sum up the individual contributions to compute value of the 
% conditional log likelihood function
CLL_value = sum(l_vector);
end