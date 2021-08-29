%  *************************************************************
%  * STEP 1: estimate VAR(2) in standard form equation-wise    *
%  *         to obtain c, Phi_1 and Phi_2                      *
%  *************************************************************
%     This procedure estimates the parameters of a VAR(2) in reduced form
%     equation-wise:
% 
%     y_t = c + Phi_1 * y_t-1 + Phi_2 * y_t-2 + epsilon_t
% 
%     where y_t, c, y_t-1, y_t-2 and epsilon_t are (n x 1) vectors
%     and Phi_1 and Phi_2 (n x n) matrices of OLS estimates
%     
%     Input arguments are:
%     x_t: Data set
%     p: number of lags to be included (1 or 2)
 
function [constant,Phi_1,Phi_2,Omega] = PHI(x_t,p)

  if p >= 3
      disp('Function is not defined for this number of lags');
      return
  else
      T = size(x_t,1);
      n = size(x_t,2);
      
      % Set up the X-matrix
      X = ones(T-p,n*p+1);        % Include a constant and the first two lags of all four system variables
      
      x_1 = lagmatrix(x_t,1);     %  w_t-1, x_t-1, y_t-1 and z_t-1 
      if p == 2
      x_2 = lagmatrix(x_t,2);     %  w_t-2, x_t-2, y_t-2 and z_t-2 
      end
      
      X(:,1:n) = x_1(p+1:size(x_t,1),:);    % first lags 
      if p == 2
      X(:,n+1:2*n) = x_2(p+1:size(x_t,1),:);% second lags 
      end

      epsilon = ones(T-p,n);                % vector for OLS residuals 
      
     %*************************************************************
     %* estimate phi parameters and the constant equation-wise    *
     %*************************************************************

      coefficients = zeros(n*p+1,n);

      for i = 1:n

      % set up y_vector: equation-wise estimation 
      y = x_t(p+1:size(x_t,1),i);     % run w_t, x_t, y_t and z_t separately 
      beta = (X'*X)\(X'*y);           % (n*p+1)x1 vector for each system variable:
                                      % for p=1:
                                      % c, phi1_1, phi2_1, phi3_1, phi4_1 
                                      % for p=2:
                                      % c, phi1_1, phi2_1, phi3_1, phi4_1, phi1_2, phi2_2, phi3_2 phi4_2  

      coefficients(:,i) = beta;
      e_t = y-X*beta;                 % residual (T-p)x1 vector of equationwise OLS
      epsilon(1:size(e_t,1),i) = e_t;
      end
      
      constant = coefficients(p*n+1,:)'; % Create vector of constant 

      Phi_1 = coefficients(1:n,:)';      % Create Phi1 (matrix of coefficients of the first lag)
      if p == 2
      Phi_2 = coefficients(n+1:2*n,:)';  % Create Phi2 (matrix of coefficients of the second lag)
      else
      Phi_2 = [];
      end

      Omega = (epsilon'*epsilon)/(T-p);   % variance-covariance matrix of composite shocks 

      
  end
      

end