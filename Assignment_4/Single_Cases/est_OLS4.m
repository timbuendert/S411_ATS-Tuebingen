function [output_vector] = est_OLS4(RW)

T = size(RW,1)-1;

X = [ones(T,1), (1:1:T)', RW(1:T)];
y = RW(2:T+1);

b_est = inv(X'*X)*(X'*y);
alpha_est = b_est(1);
delta_est = b_est(2);
rho_est = b_est(3);

s_sq = (1/(T-3))*sum((RW(2:T+1)-alpha_est-delta_est*(1:1:T)'-rho_est*(RW(1:T))).^2);
var_b = s_sq * inv(X'*X);
rho_se = sqrt(var_b(3,3));

rho_tstat = (rho_est-1)/rho_se;

rho_testst = T*(rho_est-1);

output_vector = [rho_est, rho_se, rho_testst, rho_tstat];
end