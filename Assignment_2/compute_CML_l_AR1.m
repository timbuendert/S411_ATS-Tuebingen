function [l_vector] = compute_CML_l_AR1(theta, y)
% read out parameters from theta
c = theta(1);
phi = theta(2);
sigma_2 = theta(3);

% calculate sequence vector of epsilon_t
e_t = y(2:size(y,1)) - c - phi*y(1:(size(y,1)-1));

% compute vector of conditional log likelihood contributions
l_vector = log(1/(sqrt(2*pi*sigma_2))) - ((e_t.^2)/(2*sigma_2));


% or: 

%l_vector = zeros((size(y,1)-1), 1);
%for i = 2:size(y,1)
%    e_t = y(i) - c - phi*y(i-1);
%    l_vector(i-1, 1) = log(1/(sqrt(2*pi*sigma_2)))-(e_t^2/(2*sigma_2));
%end 

end

