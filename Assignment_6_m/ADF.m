%% Function to conduct an ADF case with k = 2 for a series "y_t" based on case "est_case"

function [rho, se_rho, t_statistic] = ADF(y_t, est_case)

% set up vector for constant
constant = ones(length(y_t), 1);

% set up vector for time trend (since we lag the difference two times, 
% we need to start at -2)
T = (-2:(length(y_t)-3))';

% set up y_t-1
y_t1 = lagmatrix(y_t, 1);

% set up delta y_t-1
delta_y_t1 = lagmatrix(y_t, 1) - lagmatrix(y_t, 2);

% set up delta y_t-2
delta_y_t2 = lagmatrix(y_t, 2) - lagmatrix(y_t, 3);



% set up X matrix depending on chosen case
if est_case == "case1"
    X = [y_t1, delta_y_t1, delta_y_t2]; 
elseif est_case == "case2"
    X = [constant, y_t1, delta_y_t1, delta_y_t2];
elseif est_case == "case4"
    X = [constant, T, y_t1, delta_y_t1, delta_y_t2]; 
else
    disp("!!! Error setting up X-matrix when choosing the case !!!")
end

% remove rows with NAs
na_rows = any(isnan(X),2);
X(na_rows,:) = [];


% set up y matrix
y = y_t(4:length(y_t)); 


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




% compute standard errors
% set up s^2 depending on the case and then calculate variance and
% obtain standard error for rho
if est_case == "case1"
    s2 = (1/(length(X)-3)*sum((y - rho*X(:,1) - beta(2)*X(:,2) - beta(3)*X(:,3)).^2));
    var_beta = s2 * inv(X'*X);
    se_rho = sqrt(var_beta(1,1));
elseif est_case == "case2"
    s2 = (1/(length(X)-4)*sum((y - beta(1) - rho*X(:,2) - beta(3)*X(:,3) - beta(4)*X(:,4)).^2));
    var_beta = s2 * inv(X'*X);
    se_rho = sqrt(var_beta(2,2));
elseif est_case == "case4"
    s2 = (1/(length(X)-5)*sum((y - beta(1) - beta(2)*X(:,2) - rho*X(:,3) - beta(4)*X(:,4) - beta(5)*X(:,5)).^2));
    var_beta = s2 * inv(X'*X);
    se_rho = sqrt(var_beta(3,3));
else
    disp("!!! Error calculating std. errors when choosing the case !!!")
end


% calculate t-statistic
t_statistic = (rho - 1)/se_rho;


end