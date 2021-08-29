%  ***************************************************************
%  * PROC 3: sequences of moving average parameters psi          *                   
%  * Vector MA(infininty) representation: psi sequence           *
%  ***************************************************************
%  * This procedure:  1.  sets up the F matrix that contains Phi_1, Phi_2 
%                     2.  decomposes F into the np eigenvalues lambda
%                         and eigenvectors T
%                     3.  decides on stationarity: check if largest eigenvalue 
%                         is smaller than 1 in abs. value or if lambdas are complex 
%                         uses the modulus |R|<1
%                     4.  calculates F^s for altering s and writes out the 
%                         blocks of F^s of interest to obtain the psis
%                     5.  reads out the cumulative psis
function [F,Psi] = VAR_2(Phi_1,Phi_2,S_max,p)

  % Define number of system variables 
  n = size(Phi_1,2);
  % Set up F matrix 
  if p == 1
	x1 = Phi_1;
	F = x1;
  elseif p == 2
	x1 = [Phi_1 Phi_2];
	x2 = eye(n);
	x3 = zeros(n,n);
	x4 = [x2 x3];
	F = [x1;x4];
  end
  % vector of eigenvalues lambda (np x 1) and matrix of eigenvectors (np x np)
  [T_mat,lambda] = eig(F); %dalia: new changed order, because it gave an error in 2018
  lambda = diag(lambda);
  % decide on stationarity: |lambda|< 1 or |R|<1
  stat = ones(n*p,1);
  for j = 1:n*p       % for all n*p eigenvalues 
	a = real(lambda(j));
    b = imag(lambda(j));
	if b == 0          % real eigenvalues (redundant, just for illustrative reasons) 
		stat(j) = abs(lambda(j));
    else               % complex eigenvalues 
		R = sqrt(a^2+b^2);
		stat(j) = abs(R);
    end
  end
stationarity = max(stat);
%  ************************************************************
%  * sequence of moving average parameters:                   *
%  * blocks of F^s matrix for s=1,...,infinity                *
%  ************************************************************
% set up np x np matrix with eigenvalues on main diagonal zero else 
diagonal = eye(n*p).*lambda;
% initialize psi matrix 
Psi = zeros(n*S_max,n);
psi_sequence = zeros(S_max,n*n);        % psi_sequence=F^s_11 
j = 0;                                  % it needs to jump by n in order not to overwrite anything 
for s = 1:S_max
	F_s = T_mat*(diagonal)^s*inv(T_mat);% compute F^s to obtain the psi_s matrix: 
                                        % F^s_11 is the nxn block in the upper left corner (first lag) 
                                        % F^s_12 the nxn block in the upper right corner (second lag)  
	Psi(1+j:n+j,:) = F_s(1:n,1:n); 
    for i = 1:n
		% fill psi_sequence=F^s_11         
		psi_sequence(s,(i-1)*n+1:i*n) = F_s(i,1:n);      % psi_11, psi_12, ..., psi_1n 
														 % psi_21, psi_22, ..., psi_2n
														 %	...     ...     ...   ...
														 %	psi_n1, psi_n2, ..., psi_nn    
    end  
j = j+n;
end
% check that all psis are real numbers i.e. imag(Psi)==0 
imaginary = zeros(size(Psi,1),size(Psi,2));
for r = 1:size(Psi,1)
    for c = 1:size(Psi,2)
		if imag(Psi(r,c)) ~= 0
			imaginary(r,c) = 1;
        end
    end
end
if sum(sum(imaginary)) ~= 0
	disp('not all psis are real numbers');
end

% remark: if the imanginary part is zero but there has been a non-zero imaginary 
%            component, this still prints the error message, even if it shouldn't.        
Psi = real(Psi);
psi_sequence = real(psi_sequence);

%  ************************************************************
%  * absolute summability of psi coefficients: elementwise    *
%  ************************************************************
element_summability = sum(abs(psi_sequence'));
abs_summability = sum(sum(abs(psi_sequence)));
% needed to plot 
x_axis = 1:S_max;

end