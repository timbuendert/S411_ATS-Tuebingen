%% Function to conduct an ADF case for "k" lags for a series "y_t" based on case "est_case"

function [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(y_t, est_case, k)

%% X matrix
% set up vector for constant
constant = ones(length(y_t), 1);

% set up vector for time trend, depending on number of lags
T = (-k:(length(y_t)-(k+1)))';

% set up y_t-1 vector
y_t1 = lagmatrix(y_t, 1);


% set up matrix for lagged differences

% set up empty matrix to store the values
delta_y_matrix = zeros(length(y_t), k);

for i = 1:k
    
    % compute lagged difference
    delta_yt = lagmatrix(y_t, i) - lagmatrix(y_t, (i+1));
    
    % store in matrix
    delta_y_matrix(:,i) = delta_yt;
    
end


% set up X matrix depending on case
if est_case == "case1"
    X = [y_t1, delta_y_matrix]; 
elseif est_case == "case2"
    X = [constant, y_t1, delta_y_matrix];
elseif est_case == "case4"
    X = [constant, T, y_t1, delta_y_matrix]; 
else
    disp("!!! Error setting up X-matrix when choosing the case !!!")
end


% remove rows with NAs
na_rows = any(isnan(X),2);
X(na_rows,:) = [];

%% y matrix

% set up y matrix depending on the lag chosen
y = y_t((k+2):length(y_t)); 

%% OLS estimation

% calculate beta (vector of parameter estimates)
beta = inv(X'*X) * (X'*y);

% obtain rho
if est_case == "case1"
    rho = beta(1); 
elseif est_case == "case2"
    rho = beta(2);
elseif est_case == "case4"
    rho = beta(3);
else
    disp("!!! Error obtaining rho when choosing the case !!!")
end

%% compute standard error of rho

% compute standard errors
% set up s^2 depending on the case and then calculate variance and
% obtain standard error for rho
if est_case == "case1"
    % set up m
    m = k + 1;
    % set up s^2
    s2 = (1/(length(X)-m)*sum((y - sum((beta' .* X),2)).^2));
    % compute estimated variance matrix of beta depending on s^2
    var_beta = s2 * inv(X'*X);
    % read out the standard error of rho
    se_rho = sqrt(var_beta(1,1));
elseif est_case == "case2"
    m = k + 2;
    s2 = (1/(length(X)-m)*sum((y - sum((beta' .* X),2)).^2));
    var_beta = s2 * inv(X'*X);
    se_rho = sqrt(var_beta(2,2));
elseif est_case == "case4"
    m = k + 3;
    s2 = (1/(length(X)-m)*sum((y - sum((beta' .* X),2)).^2));
    var_beta = s2 * inv(X'*X);
    se_rho = sqrt(var_beta(3,3));
else
    disp("!!! Error calculating std. errors when choosing the case !!!")
end

%% calculate t-statistic and SBC / AIC
% calculate t-statistic
t_statistic = (rho - 1)/se_rho;

% calculate SBC and AIC
SBC = log(s2) + (m/length(X))*log(length(X));
AIC = log(s2) + (2*m)/length(X);


end