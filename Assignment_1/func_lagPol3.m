function [evaluations] = func_lagPol3(phi1, phi2, phi3, input_val)
% function to compute evaluations of given lag polynomial for input values

evaluations = zeros(size(input_val, 2), 1); % initialize evaluations vector
for i = 1:size(input_val, 2) % iterate through all input values
    z = input_val(i); % determine respective input value
    
    % evaluate lap polynomial based on z
    evaluations(i,1) = 1 - (phi1*z) - (phi2*(z)^2) - (phi3*(z)^3);
end
end

