function [output_vector] = est_OLS(RW, case_t)

T = size(RW,1)-1; % determine length of realization

% determine X based on case
if case_t == 1
    X = RW(1:T);
elseif case_t == 2
    X = [ones(T,1), RW(1:T)];
elseif case_t == 4
    X = [ones(T,1), (1:1:T)', RW(1:T)];
end

y = RW(2:T+1); % define y

b_est = inv(X'*X)*(X'*y); % estimate RW parameters

% retrieve parameter estimates based on case
if case_t == 1
    rho_est = b_est(1);
elseif case_t == 2
    alpha_est = b_est(1);
    rho_est = b_est(2);    
elseif case_t == 4
    alpha_est = b_est(1);
    delta_est = b_est(2);
    rho_est = b_est(3);   
end

% determine se_rho based on case
if case_t == 1
    s_sq = (1/(T-1))*sum((RW(2:T+1)-rho_est*RW(1:T)).^2);
    var_rho = s_sq * inv(X'*X);
    rho_se = sqrt(var_rho);
elseif case_t == 2
    s_sq = (1/(T-2))*sum((RW(2:T+1)-alpha_est-rho_est*(RW(1:T))).^2);
    var_b = s_sq * inv(X'*X);
    rho_se = sqrt(var_b(2,2));
elseif case_t == 4
    s_sq = (1/(T-3))*sum((RW(2:T+1)-alpha_est-delta_est*(1:1:T)'-rho_est*(RW(1:T))).^2);
    var_b = s_sq * inv(X'*X);
    rho_se = sqrt(var_b(3,3));
end

% calculate t statistic
rho_tstat = (rho_est-1)/rho_se;

% calculate test statistic: T * (ρ_est - 1)
rho_testst = T*(rho_est-1);

% collect objects in output vector
output_vector = [rho_est, rho_se, rho_testst, rho_tstat];
end