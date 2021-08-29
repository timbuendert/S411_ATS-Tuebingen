function [output_vector] = est_OLS1(RW)

T = size(RW,1)-1;

X = RW(1:T);
y = RW(2:T+1);

rho_est = inv(X'*X)*(X'*y);

s_sq = (1/(T-1))*sum((RW(2:T+1)-rho_est*(RW(1:T))).^2);
var_rho = s_sq * inv(X'*X);
rho_se = sqrt(var_rho);

rho_tstat = (rho_est-1)/rho_se;

rho_testst = T*(rho_est-1);

output_vector = [rho_est, rho_se, rho_testst, rho_tstat];
end

