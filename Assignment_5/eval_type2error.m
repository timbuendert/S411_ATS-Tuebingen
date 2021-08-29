function [test_error,t_error] = eval_type2error(test_RW, t_RW, test_st, t_st)
% type II error: if not reject the incorrect null hypothesis of unit root

% consider one-sided test with α = 0.05: compute critical values
% (i.e. 5%-quantiles)
test_crit = quantile(test_RW, 0.05);
t_crit = quantile(t_RW, 0.05);

% count errors: where statistic > critical value -> not reject false H0
% divide by length to compute relative frequency
test_error = sum(test_st > test_crit)/length(test_st);
t_error = sum(t_st > t_crit)/length(t_st);

end