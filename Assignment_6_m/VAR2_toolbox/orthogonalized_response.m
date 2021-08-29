%  ***************************************************************
%  * PROC 4: orthogonalized impulse response functions:          *
%  *         response to pure innovations                        *
%  *         plot sequences of moving average parameters psi*aj  *                   
%  ***************************************************************
%     This procedure calculates a weighted Psi-matrix sequence 
%     that can be used to plot orthogonalized impulse response functions:
%             - we use a Cholesky decomposition to back out the structure: A and D
%             - sd=0: implemented for one unit shocks in pure innovations
%                     weighting matrix used: weight = A
%             - sd=1: implemented for one standard deviation shocks in u_t
%                     weighting matric used: weight = P=A*C=A* sqrt(D)
%             - we calculate the modified Psi-sequence by: Psi_s * weight         
%             - psi_orthogonalizied is a matrix with all elements of Psi_s for all s
%               as columns (n*n columns) and s rows
% 	Input arguments: 
% 			Psi: taken from procedure VAR_2
% 			A and C: taken from procedure Cholesky decomposition
% 			S_max: maximum number of lags considered
% 			sd: 0 if one unit shocks are considered
% 				1 if one standard deviation shocks are considered
% 	Output arguments:
% 			matrix of orthogonalized impulse responses.
% 			The first column refers to the effect of a shock in the first variable on
% 			the first variable. The second one refers to a shock in the second variable 
% 			on the first variable, and so on. The (n+1)th column refers to the effect
% 			of a shock in the first variable on the second one...
function psi_orthogonalized = orthogonalized_response(Psi,A,C,S_max,sd)

% Define number of system variables 
n = size(Psi,2);
Psi_s = ones(n,n);
psi_orthogonalized = zeros(S_max+1,n*n);             % psi_ortho=Psi*A               
j = 0; % it needs to jump by n in order not to overwrite anything 
% into first row s=1  and sd=1: P = A * C = A * sqrt(D)  
%     for variable 1: sd(u_1t)        , 0               , 0               , 0
%     for variable 2: a_21 * sd(u_1t) , sd(u_2t)        , 0               , 0 
%     for variable 3: a_31 * sd(u_1t) , a_32 * sd(u_2t) , sd(u_3t)        , 0
%     for variable 4: a_41 * sd(u_1t) , a_42 * sd(u_2t) , a_43 * sd(u_3t) , sd(u_4t) 
%     into first row s=1  and sd=0: A                                                  
if sd == 0               % one unit shock in pure innovation                        
    Psi_0 = A;
    disp('one unit shock in pure innovations');
elseif sd == 1          % one standard deviation shock in pure innovations         
    Psi_0 = A*C;
    disp('one standard deviation shock in pure innovations');
end
for i = 1:n
	psi_orthogonalized(1,(i-1)*n+1:i*n) = Psi_0(i,1:n);         
end
for s = 1:S_max
	Psi_s = Psi(1+j:n+j,:);
	if sd == 0
		orthogonal_psi = Psi_s*A;
	elseif sd == 1 
		orthogonal_psi = Psi_s*A*C;
    end
	for i = 1:n
		% fill psi*_sequence=F^s_11         
		psi_orthogonalized(s+1,(i-1)*n+1:i*n)=orthogonal_psi(i,1:n);         
		% 	psi*_11, psi*_12, ..., psi*_1n 
		%	psi*_21, psi*_22, ..., psi*_2n
		%	...     ...     ...   ...
		%	psi*_n1, psi*_n2, ..., psi*_nn    
    end     
	j = j+n;
end


end