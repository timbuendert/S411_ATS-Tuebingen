function [L_value] = CLL_value(theta, y)
l = compute_CML_l_AR1(theta, y); % compute log likelihood contributions vector based on parameters and data
L_value = sum(l); % sums l up to get value of CLL
end

