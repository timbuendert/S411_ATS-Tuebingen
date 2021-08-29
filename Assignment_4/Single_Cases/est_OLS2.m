function [output_vector] = est_OLS2(RW)

T = size(RW,1)-1;

X = [ones(T,1), RW(1:T)];
y = RW(2:T+1);

b_est = inv(X'*X)*(X'*y);
alpha_est = b_est(1);
rho_est = b_est(2);

s_sq = (1/(T-2))*sum((RW(2:T+1)-alpha_est-rho_est*(RW(1:T))).^2);
var_b = s_sq * inv(X'*X);
rho_se = sqrt(var_b(2,2));

rho_tstat = (rho_est-1)/rho_se;

rho_testst = T*(rho_est-1);

output_vector = [rho_est, rho_se, rho_testst, rho_tstat];
end