function [l_vector] = lik_contrib(param, y)
% lik_contrib(param, y): function to compute the conditional log-likelihood
% contributions based on a column vector of parameters (param) and data
% (y). It returns a (Tx1) column vector of the individual log likelihood
% contributions (l_vector).

% read out parameters from parameter vector
mu = param(1);
theta = param(2);
sigma_2 = param(3);


% calculate sequence vector of ɛ_t in an iterative approach
e_t = zeros(size(y,1), 1); % define vector of size (Tx1) storing the contributions

for i = 1:size(y,1) % iterate from 1 to T
    if i == 1
        e_t(i,1) = y(i) - mu - theta * 0; % assuming ɛ_0 = 0
    else
        e_t(i,1) = y(i) - mu - theta * e_t(i-1); % determine ɛ_t
    end
end


% compute vector of conditional log likelihood contributions by
% vectorization (element-wise squaring the sequence vector of ɛ_t)
l_vector = log(1/(sqrt(2*pi*sigma_2))) - ((e_t.^2)/(2*sigma_2));

end