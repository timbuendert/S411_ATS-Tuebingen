%  ***************************************************************
%  * PROC 5: MSE or variance decomposition Omega= A D A'         *                   
%  ***************************************************************
%     This procedure performs a variance decomposition as illustrated on page 263 (268) in the script.
% 	Input arguments: 
% 			Psi: taken from procedure VAR_2
% 			A and D: taken from procedure Cholesky decomposition
% 			S_max: maximum number of lags considered
% 	Output arguments:
% 			matrix of variance decompositions
% 			The first column refers to the fraction of the variance of variable 1 that is due to 1.
% 			The second column refers to the fraction of the variance of variable 2 that is due to 1.
% 			The (n+1)th column refers to the fraction of the variance of variable 1 that is due to 2, etc. 
function variance_dec = variance_deco(Psi,A,D,S_max)
% Define number of system variables
n = size(Psi,2);
% Initialize some matrices
MSE_mat = zeros(n*S_max,n*n); 
variance_dec = zeros(S_max,n*n);
MSE_diag = zeros(S_max,n*n);
Var_u_t = diag(D);        % read out variances of idiosyncratic shocks 
%  The matrix MSE_mat is of dimensions (S_max*n)x(n*n). The first nxn block refers to the 
%  contributions to the MSE of the 1-period-ahead forecast that are related to the first 
%  variable. The [1:n,n+1:2n] elements of MSE_mat represent contributions to the MSE of the
%  1-step-ahead forecast that are related to the second variable etc.  
for s = 1:S_max
% compute matrix of MSE using matrices Omega and Psi_1 to Psi_S-1 
    if s == 1 % compute MSE contributions for 1-step-ahead forecast
		for j = 1:n
			MSE_mat(1:n,n*(j-1)+1:n*j) = Var_u_t(j)*A(:,j)*A(:,j)'; % see page 263 in the script
        end
    else % compute MSE contributions for more steps ahead
		Psi_mat = Psi(n*(s-2)+1:n*(s-1),:);  % Read out respective Psi_s
		for j = 1:n  % see page 263 in the script
			MSE_mat(n*(s-1)+1:n*s,n*(j-1)+1:n*j) = MSE_mat(n*(s-2)+1:n*(s-1),n*(j-1)+1:n*j)+Var_u_t(j)*Psi_mat*A(:,j)*A(:,j)'*Psi_mat';
        end
    end
end
for j = 1:S_max
for i = 1:n
% read out diagonal elements of the (nxn)-size blocks in MSE_mat
MSE_diag(j,n*(i-1)+1:n*i) = diag(MSE_mat(n*(j-1)+1:n*j,n*(i-1)+1:n*i))';	
end
end
for i = 1:n
% 	 Compute elements on the main diagonal of the MSE matrix (first row of MSE_sum
%    refers to s=1, second row to s=2,...
	if i == 1
	MSE_sum = MSE_diag(:,n*(i-1)+1:n*i); 
    else
	MSE_sum=MSE_sum+MSE_diag(:,n*(i-1)+1:n*i);
    end
end
for i = 1:n
% 	Compute the fraction that each of the single variables contributed to the elements on
% 	the main diagonal of the MSE matrix (for each s) 
	variance_dec(:,n*(i-1)+1:n*i) = MSE_diag(:,n*(i-1)+1:n*i)./MSE_sum;
end
end