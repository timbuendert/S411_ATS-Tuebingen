function [L_value] = CLL_value_min(theta, y)
l = compute_CML_l_AR1(theta, y); % compute log likelihood contributions vector based on parameters and data
L_value = -(sum(l)); % sum l up to get value of CLL and take opposite sign -> as searching for minimum with functions
end