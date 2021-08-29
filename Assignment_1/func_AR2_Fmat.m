function [Fmat, Feig] = func_AR2_Fmat(phi1,phi2,phi3)
% function to create a 3x3 F matrix: given phi coefficients in  first row, 
% set last column all 0 and remaining part as identity matrix

Fmat = [phi1 phi2 phi3; eye(2,3)];
Feig = eig(Fmat); % compute the eigenvalues of F
end