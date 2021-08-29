% Assignment #3 by 
% Leonard Berger (Student-ID: 5629475) 
% Tim-Moritz Bündert (Student-ID: 5635975)

% The following function files accompany this code: 
% - lik_contrib.m
% - lik_func.m
% - lik_func_min.m  
% - func_MA1.m
% They are explained in more detail in their respective files and below 
% when they are used.

clear
clc
close all

% Regarding the options used in the parameter estimation with the CML 
% toolbox: to produce a cleaner command window, the first two arguments 
% 'Display' & 'iter' can be omitted (see option 1). Alternatively, option 2
% presents the orginal options given in the Assignment task.

% 1) options = optimset('TolX',10e-40,'TolFun',10e-40, 'MaxIter',10^10 , 'MaxFunEvals', 100000);
% 2) options = optimset('Display','iter','TolX',10e-40,'TolFun',10e-40, 'MaxIter',10^10 , 'MaxFunEvals', 100000);

% In case option 1 is preferred, please specify "opt_est = 1;"
% In case option 2 is preferred, please specify "opt_est = 2;"
opt_est = 1;

if opt_est == 1
    options = optimset('TolX',10e-40,'TolFun',10e-40, 'MaxIter',10^10 , 'MaxFunEvals', 100000);
else
    options = optimset('Display','iter','TolX',10e-40,'TolFun',10e-40, 'MaxIter',10^10 , 'MaxFunEvals', 100000);
end


rng(1402); % set seed to ensure replicability 

%% Exercise 1-1 (Section 2 in paper)
T = 100; % set number of observations

% func_MA1(μ, θ, σ^2, T): function to generate a realization of an MA(1)
% process based on the supplied parameters. It returns a vector of size T
% containing the respective process values.

% simulate MA(1) process based on parameters: μ=10, θ=0.2 σ^2=1
MA1_realiz = func_MA1(10, 0.2, 1, T); 

% plot the simulated MA(1) realization
figure('Position', [0 0 650 450]); % fix size of plot
plot(1:T, MA1_realiz, 'k', "linewidth", 1); % plot realization of process
yline(10, '--r', "linewidth", 1) % add E[MA(1)]
legend("MA(1) realization", "$\mu$","location", "northeast",'Interpreter','latex', 'Fontsize', 13);
title('Realization of an MA(1) process');
xlabel({'','$t$'}, 'Fontsize', 15, 'Interpreter',"latex"); % create whitespace between axis and label by multiline label
ylabel({'','$y_{t}$'}, 'Fontsize', 15, 'Interpreter',"latex"); % create whitespace between axis and label by multiline label
axis([0 T 7 14]);
%saveas(gcf,'MA1Realization.png') % optionally: save plot


%% Exercise 1-2 (Section 3.1 in paper)
% lik_contrib(param, y): function to compute the individual conditional 
% log-likelihood contributions based on a column vector of parameters 
% (param) and data (y). It returns a (Tx1) column vector of the individual 
% log likelihood contributions.

%% Exercise 1-3 (Section 3.1 in paper)
fprintf('Exercise 1-3 (Section 3.1 in paper) ------------------------------------------------');

% lik_func(param, y): function to first, determine the vector of individual
% log-likelihood contributions based on parameters (param) and data
% (y), and, second, summing up all individual contributions to compute and 
% return the value of the total conditional log-likelihood function
% (CLL_value)

% compute value of the conditional log likelihood (CLL) function for 
% parameter specifications (a) and (b)
CLL_a = lik_func([10, 0.2, 1]', MA1_realiz);
CLL_b = lik_func([5, 0.7, 2]', MA1_realiz);

fprintf('\n\nCLL value for parameters (a): %.3f\nCLL value for parameters (b): %.3f\n',[CLL_a, CLL_b]')

% Interpretation:
fprintf('\n- CLL with parameters of (a) results in higher CLL value\n')
fprintf('  -> reasonable as parameters of (a) are the true parameters of data-generating process')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 2-1
addpath('CML'); % include the CML toolbox for the analysis

%% Exercise 2-2 (Section 3.2 in paper)
% Overview of different algorithms and estimates for the
% variance-covariance matrix in the CML toolbox:
% algorithm:    1 - fminsearch
%               2 - fminunc
%               3 - patternsearch
% covPar:       1 - Hessian-based covariance matrix
%               2 - OPG-based covariance matrix
%               3 - QML-based covariance matrix

% lik_func_min(param, y): same as @lik_func, but returns the negative value 
% of the total conditional log-likelihood function because the optimization
% algorithms are minimization algorithms

% Using @lik_func_min, @lik_contrib to estimate μ, θ, and σ^2 for the MA(1)
% realization (MA1_realiz) by conditional ML. The starting values are the
% true values and the algorithm 1 (fminsearch) and covariance matrix 1 
% (Hessian-based) are selected (with the defined options).
[param_est_11, CLL_value_11, g_11, cov_11, retcode_11] = CML(@lik_func_min, @lik_contrib, MA1_realiz, [10, 0.2, 1]', 1, 1, options);


%% Exercise 2-3 (Section 3.3 in paper)
fprintf('\n\n');
fprintf('Exercise 2-3 (Section 3.3 in paper) ------------------------------------------------');

% The standard errors are computed by taking the square root of the kth
% row, kth column element of the estimated covariance matrix
se_mu_11 = sqrt(cov_11(1,1));
se_theta_11 = sqrt(cov_11(2,2));
se_sigma2_11 = sqrt(cov_11(3,3));

% Display the standard errors together with the true and estimated
% parameter values
fprintf('\n\n');
fprintf('                          μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f \n',[["True              ";"Estimate (Ex. 2-2)"; "Standard Error    "], [10; param_est_11(1); se_mu_11], [0.2; param_est_11(2); se_theta_11], [1; param_est_11(3); se_sigma2_11]]')

%% Exercise 2-4 (Section 3.4 in paper)
fprintf('\n\n');
fprintf('Exercise 2-4 (Section 3.4 in paper) ------------------------------------------------');

% Consider H_0: θ = 0.4: compute t-statistic by (θ_hat - θ_bar) / s.e.(θ_hat)
t_theta_11 = (param_est_11(2) - 0.4)/ se_theta_11;

% compute critical value of SND given two-sided test with α = 0.05 (because
% t ~ N(0,1))
x_025_SND = norminv(0.05/2);

% Interpretation of result
fprintf('\n\nt = %.3f\n\n',t_theta_11)
fprintf('0.05/2 = 0.025 quantile of SND = %.3f\n', x_025_SND)
fprintf('1-(0.05/2) = 0.975 quantile of SND = %.3f\n\n', -x_025_SND)

fprintf('As %.3f < t = %.3f < %.3f:\n', [x_025_SND, t_theta_11, -x_025_SND])
fprintf('H_0: θ = 0.4 cannot be rejected at α = 0.05 \n')

%% Exercise 2-5 (Section 3.5 in paper)
fprintf('\n\n');
fprintf('Exercise 2-5 (Section 3.5 in paper) ------------------------------------------------');

% compute probability to the right of |t-statistic| and multiply by two
% to obtain two-sided p-value of t from Exercise 2-4
p_theta_11 = 2*(1-normcdf(abs(t_theta_11)));

% Interpretation of result
fprintf('\n\nThe two-sided p-value of t (from 2-4) is %.3f.\n\n', p_theta_11) 
fprintf('As %.2f > %.2f = α, the H_0: θ = 0.4 cannot be rejected at α = 0.05.\n', p_theta_11, 0.05)

%% Exercise 2-6 (Section 3.6 in paper)
fprintf('\n\n');
fprintf('Exercise 2-6 (Section 3.6 in paper) ------------------------------------------------');

% Computation of 95% confidence interval by:
% lower/upper bound = μ_hat +/- |0.025 quantile of SND| * s.e.(μ_hat)
c_low_mu_11 = param_est_11(1) - abs(x_025_SND) * se_mu_11;
c_upp_mu_11 = param_est_11(1) + abs(x_025_SND) * se_mu_11;

fprintf('\n\nμ:')
fprintf('\nThe Confidence Interval ranges from %.3f to %.3f (estimate: %.3f). \n',[c_low_mu_11, c_upp_mu_11, param_est_11(1)])
fprintf('1) For any μ_bar in [%.3f, %.3f], the H_0: μ_0 = μ_bar cannot be rejected at α = 0.05. \n',[c_low_mu_11, c_upp_mu_11])
disp('2) When analysing different realizations of the process and its parameters, in 95% of the cases (on average) μ_0 lies within the bounds of the stoachstic confidence interval.')

% Computation of 95% confidence interval by:
% lower/upper bound = θ_hat +/- |0.025 quantile of SND| * s.e.(θ_hat)
c_low_theta_11 = param_est_11(2) - abs(x_025_SND) * se_theta_11;
c_upp_theta_11 = param_est_11(2) + abs(x_025_SND) * se_theta_11;

fprintf('\nθ:')
fprintf('\nThe Confidence Interval ranges from %.3f to %.3f (estimate: %.3f). \n',[c_low_theta_11, c_upp_theta_11, param_est_11(2)])
fprintf('1) For any θ_bar in [%.3f, %.3f], the H_0: θ_0 = θ_bar cannot be rejected at α = 0.05. \n',[c_low_theta_11, c_upp_theta_11])
disp('2) When analysing different realizations of the process and its parameters, in 95% of the cases (on average) θ_0 lies within the bounds of the stoachstic confidence interval.')

fprintf('\nHowever, it is important to mention that interpretation 2) is not very useful in this regard, because only a single realization is considered.\n')

% Add-On to Exercise 2-6:
fprintf('\nIn the present case, θ_0 is not captured by the CI.\n')
fprintf('This might be due to Interpretation 2) of the confidence intervals.\n')
fprintf('To evaluate this statement, 1000 estimations and CI are computed:\n')


n_iter = 1000; 

% define 1000 (n_iter) different random numbers which are used as seeds
seeds = randsample(n_iter*50,n_iter,false)';

range = (1:1:n_iter)'; % define range for iterating

% initiate empty vectors
c_low = zeros(n_iter,1);
c_upp = zeros(n_iter,1);
est_t = zeros(n_iter,1);
decision = zeros(n_iter,1);

ind = 1; % initiate index

for i = seeds
    rng(i); % set seed to ensure replicability
    MA1_r = func_MA1(10, 0.2, 1, 100); % simulate MA(1) process based on parameters
    
    % estimation of parameters
    [param_est, ~, ~, cov, ~] = CML(@lik_func_min, @lik_contrib, MA1_r, [10, 0.2, 1]', 1, 1, options);

    % construction of the confidence interval
    c_low_t = param_est(2) - abs(x_025_SND) * sqrt(cov_11(2,2));
    c_upp_t = param_est(2) + abs(x_025_SND) * sqrt(cov_11(2,2));

    % evaluate whether 0.2 (θ_0) lies in the confidence interval
    if 0.2 >= c_low_t  && 0.2 <= c_upp_t
        dec_eval = 1;
    else
        dec_eval = 0;
    end 
    
    % save results of this iteration in the respective vector
    c_low(ind, 1) = c_low_t;
    c_upp(ind, 1) = c_upp_t;
    est_t(ind, 1) = param_est(2);
    decision(ind, 1) = dec_eval;
    
    ind = ind + 1;
end

% optionally: print for each iteration the lower and upper bound of the CI
%fprintf('\n\nθ_0 = 0.2\n\n')
%fprintf('#    C_low    C_upp     θ_0 in CI?\n');
%fprintf('%i %8.2f %8.2f %10.i\n', [range, c_low, c_upp, decision]')

fprintf('\nIn total, θ_hat is captured by the confidence interval in %i of %i iterations (%.2f percent) .\n', sum(decision), n_iter, (100*sum(decision)/n_iter));


%% Exercise 2-7 (Section 3.7 in paper)
fprintf('\n\n');
fprintf('Exercise 2-7 (Section 3.7 in paper) ------------------------------------------------');

% estimation of parameters with OPG-based variance-covariance matrix
[param_est_12, CLL_value_12, g_12, cov_12, retcode_12] = CML(@lik_func_min, @lik_contrib, MA1_realiz, [10, 0.2, 1]', 1, 2, options);
% compute respective standard errors
se_mu_12 = sqrt(cov_12(1,1));
se_theta_12 = sqrt(cov_12(2,2));
se_sigma2_12 = sqrt(cov_12(3,3));

% estimation of parameters with QML-based variance-covariance matrix
[param_est_13, CLL_value_13, g_13, cov_13, retcode_13] = CML(@lik_func_min, @lik_contrib, MA1_realiz, [10, 0.2, 1]', 1, 3, options);
% compute respective standard errors
se_mu_13 = sqrt(cov_13(1,1));
se_theta_13 = sqrt(cov_13(2,2));
se_sigma2_13 = sqrt(cov_13(3,3));

% Display standard erros with all three estimates for variance-covariance
% matrix 
fprintf('\n\nStandard Errors:\n\n');
fprintf('                     μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f \n',[["Hessian-based";"OPG-based    ";"QML-based    "], [se_mu_11; se_mu_12; se_mu_13], [se_theta_11; se_theta_12; se_theta_13], [se_sigma2_11; se_sigma2_12; se_sigma2_13]]')

% Interpretation:
fprintf('\n- close standard errors for all three alternatives\n');
fprintf('  -> indicates that Information-Matrix equality holds + conditional density is correctly specified (as it is the case)\n');


%% Exercise 2-8 (Section 3.8 in paper)
fprintf('\n\n');
fprintf('Exercise 2-8 (Section 3.8 in paper) ------------------------------------------------');

% Repeat the previous computations based on a simulated process with 
% T = 50,000

rng(1402); % set seed to ensure replicability

T_ls = 50000; % set number of observations
MA1_realiz_ls = func_MA1(10, 0.2, 1, T_ls); % simulate MA(1) process based on parameters


% Part 2

% Use @lik_func_min, @lik_contrib to estimate μ, θ, and σ^2 for the MA(1)
% realization (MA1_realiz) by conditional ML. The starting values are the
% true values and the algorithm 1 (fminsearch) and covariance matrix 1 
% (Hessian-based) are selected (with the defined options).
[param_est_11_ls, CLL_value_11_ls, g_11_ls, cov_11_ls, retcode_11_ls] = CML(@lik_func_min, @lik_contrib, MA1_realiz_ls, [10, 0.2, 1]', 1, 1, options);


% Part 3
fprintf('\n\n');
fprintf('Part 3 ------------------------------------------------');

% The standard errors are computed by taking the square root of the kth
% row, kth column element of the estimated covariance matrix
se_mu_11_ls = sqrt(cov_11_ls(1,1));
se_theta_11_ls = sqrt(cov_11_ls(2,2));
se_sigma2_11_ls = sqrt(cov_11_ls(3,3));

% Display the standard errors together with the estimated parameter values
fprintf('\n\n');
fprintf('                          μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f \n',[["True              ";"Estimate (Ex. 2-2)"; "Standard Error    "], [10; param_est_11_ls(1); se_mu_11_ls], [0.2; param_est_11_ls(2); se_theta_11_ls], [1; param_est_11_ls(3); se_sigma2_11_ls]]')


% Part 4
fprintf('\n\n');
fprintf('Part 4 ------------------------------------------------');

% Consider H_0: θ = 0.4: compute t-statistic by (θ_hat - θ_bar) / s.e.(θ_hat)
t_theta_11_ls = (param_est_11_ls(2) - 0.4)/ se_theta_11_ls;

% compute critical value of SND given two-sided test with α = 0.05 (because
% t ~ N(0,1))
x_025_SND_ls = norminv(0.05/2);

% Interpretation of result
fprintf('\n\nt = %.3f\n\n',t_theta_11_ls)
fprintf('0.05/2 = 0.025 quantile of SND = %.3f\n', x_025_SND_ls)
fprintf('1-(0.05/2) = 0.975 quantile of SND = %.3f\n\n', -x_025_SND_ls)

fprintf('As t = %.3f < %.3f < %.3f:\n', [t_theta_11_ls, x_025_SND_ls, -x_025_SND_ls])
fprintf('H_0: θ = 0.4 can be rejected at α = 0.05 \n')


% Part 5
fprintf('\n\n');
fprintf('Part 5 ------------------------------------------------');

% compute probability to the right of |t-statistic| and multiply by two
% to obtain two-sided p-value of t from exercise 2-4
p_theta_11_ls = 2*(1-normcdf(abs(t_theta_11_ls)));

fprintf('\n\nThe two-sided p-value of t (from 2-4) is %.3f.\n\n', p_theta_11_ls) 
fprintf('As two-sided p_value = %.2f < %.2f = α, the H_0: θ = 0.4 can be rejected at α = 0.05.\n', p_theta_11_ls, 0.05)


% Part 6
fprintf('\n\n');
fprintf('Part 6 ------------------------------------------------');

% Computation of 95% confidence interval by:
% lower/upper bound = μ_hat +/- |0.025 quantile of SND| * s.e.(μ_hat)
c_low_mu_11_ls = param_est_11_ls(1) - abs(x_025_SND_ls) * se_mu_11_ls;
c_upp_mu_11_ls = param_est_11_ls(1) + abs(x_025_SND_ls) * se_mu_11_ls;

% Computation of 95% confidence interval by:
% lower/upper bound = θ_hat +/- |0.025 quantile of SND| * s.e.(θ_hat)
c_low_theta_11_ls = param_est_11_ls(2) - abs(x_025_SND_ls) * se_theta_11_ls;
c_upp_theta_11_ls = param_est_11_ls(2) + abs(x_025_SND_ls) * se_theta_11_ls;

% Interpretation
fprintf('\n\nμ:')
fprintf('\nThe Confidence Interval ranges from %.3f to %.3f (estimate: %.3f). \n',[c_low_mu_11_ls, c_upp_mu_11_ls, param_est_11_ls(1)])
fprintf('1) For μ_bar in [%.3f, %.3f], the H_0: μ_0 = μ_bar cannot be rejected at α = 0.05. \n',[c_low_mu_11_ls, c_upp_mu_11_ls])
disp('2) When analysing different realizations of the process and its parameters, in 95% of the cases (on average) μ_0 lies within the bounds of the stoachstic confidence interval.')

fprintf('\n\nθ:')
fprintf('\nThe Confidence Interval ranges from %.3f to %.3f (estimate: %.3f). \n',[c_low_theta_11_ls, c_upp_theta_11_ls, param_est_11_ls(2)])
fprintf('1) For θ_bar in [%.3f, %.3f], the H_0: θ_0 = θ_bar cannot be rejected at α = 0.05. \n',[c_low_theta_11_ls, c_upp_theta_11_ls])
disp('2) When analysing different realizations of the process and its parameters, in 95% of the cases (on average) θ_0 lies within the bounds of the stoachstic confidence interval.')

fprintf('\nHowever, it is important to mention that interpretation 2) is not very useful in this regard, because only a single realization is considered.')

% Part 7
fprintf('\n\n');
fprintf('Part 7 ------------------------------------------------');

% estimation of parameters with OPG-based variance-covariance matrix
[param_est_12_ls, CLL_value_12_ls, g_12_ls, cov_12_ls, retcode_12_ls] = CML(@lik_func_min, @lik_contrib, MA1_realiz_ls, [10, 0.2, 1]', 1, 2, options);
% compute respective standard errors
se_mu_12_ls = sqrt(cov_12_ls(1,1));
se_theta_12_ls = sqrt(cov_12_ls(2,2));
se_sigma2_12_ls = sqrt(cov_12_ls(3,3));

% estimation of parameters with QML-based variance-covariance matrix
[param_est_13_ls, CLL_value_13_ls, g_13_ls, cov_13_ls, retcode_13_ls] = CML(@lik_func_min, @lik_contrib, MA1_realiz_ls, [10, 0.2, 1]', 1, 3, options);
% compute respective standard errors
se_mu_13_ls = sqrt(cov_13_ls(1,1));
se_theta_13_ls = sqrt(cov_13_ls(2,2));
se_sigma2_13_ls = sqrt(cov_13_ls(3,3));

fprintf('\n\nStandard Errors:\n\n');
fprintf('                     μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f \n',[["Hessian-based";"OPG-based    ";"QML-based    "], [se_mu_11_ls; se_mu_12_ls; se_mu_13_ls], [se_theta_11_ls; se_theta_12_ls; se_theta_13_ls], [se_sigma2_11_ls; se_sigma2_12_ls; se_sigma2_13_ls]]')

% Interpretation:
fprintf('\n- same rounded standard erros for all three alternatives')
fprintf('\n -> indicates that Information-Matrix equality holds')


fprintf('\n\nSummary of interpretation\n');

fprintf('\nWith T = 50,000 instead of T = 100:\n');
fprintf('- Estimates get very close to true parameters, while s.e. get close to 0\n');
fprintf('- H_0: θ_0 = 0.4 can be rejected at α = 0.05 with t = -46.75 -> two-sided p-value = 0\n');
fprintf('- CI for μ: [9.989, 10.010] // CI for θ: [0.187, 0.204]\n');
fprintf('- s.e. are the same for all alternatives of variance-covariance matrices across all parameters -> Information-Matrix equality holds \n');

%% Exericse 2-9 (Section 3.9 in paper)
fprintf('\n\n');
fprintf('Exercise 2-9 (Section 3.9 in paper) ------------------------------------------------');

rng(1402); % set seed to ensure replicability

% simulate MA(1) process based on given parameters
MA1_realiz_ex9 = func_MA1(10, 0.2, 1, 200); 

% set range for theta: from θ = −0.6, increase θ by 0.05 until θ = 0.8
theta_range = -0.6:0.05:0.8;

% initiate vector for CLL values of varying θ
theta_CLLval = zeros(size(theta_range, 1),1); 

ind = 1; % initiate index
for i = theta_range
    
    % Fix μ = 10 and σ^2 = 1 and compute value of the conditional log-
    % likelihood function
    L_phi = lik_func([10, i, 1], MA1_realiz_ex9);
    
    theta_CLLval(ind, 1) = L_phi; % store CLL value in vector
    ind = ind + 1; % update index
end 

% estimate parameters of MA(1) with fminsearch and Hessian-based Cov
[param_est_11_ex9,~,~,cov_11_ex9] = CML(@lik_func_min, @lik_contrib, MA1_realiz_ex9, [10, 0.2, 1]', 1, 1, options);

% plot the conditional log likelihood function for theta_range + θ_hat
% + θ_0
figure('Position', [0 0 650 450]); % fix size of plot
plot(theta_range, theta_CLLval, 'k', "linewidth", 1); % plot CLL function of process
xline(param_est_11_ex9(2), '--r', "linewidth", 1) % add CML estimate of θ (θ_hat)
xline(0.2, '--g', "linewidth", 1) % add true value of θ (θ_0)
title('Conditional log-likelihood function');
xlabel({'','$\theta$'}, 'Fontsize', 15, 'Interpreter',"latex"); % create whitespace between axis and label by multiline label
ylabel({'','$\ln \mathcal{L}$'}, 'Fontsize', 15, 'Interpreter',"latex"); % create whitespace between axis and label by multiline label
legend("$\ln \mathcal{L}(\theta)$", "$\hat{\theta}$", "$\theta_0$","location", "southeast",'Interpreter','latex', 'Fontsize', 13);
xlim([-0.6 0.8]);
%saveas(gcf,'CLLfunction_Ex9.png') % optionally: save plot


% Interpretation of Graph:
fprintf('\nSee Figure 2:');
fprintf('\n- with μ = 10 and σ^2 = 1 (true parameters): maximum CLL value at θ = 0.3207 (value of CLL estimate)\n');
fprintf('  -> small deviation from θ_0 = 0.2 \n');
fprintf('- rather peaked log likelihood function (considering y-axis) -> high estimation precision\n');
fprintf('  (as also indicated by rather low s.e.(θ_hat) = %.3f)\n', sqrt(cov_11_ex9(2,2)))


%% Exericse 2-10 (Section 3.10 in paper)
fprintf('\n\n');
fprintf('Exercise 2-10 (Section 3.10 in paper) ------------------------------------------------');

% simulate MA(1) process based on given parameters
MA1_realiz_ex10 = func_MA1(10, 5, 1, 200); 

% estimate parameters of MA(1) using fminunc and selecting the Hessian-
% based variance-covariance estimation
% select true values of μ and σ^2 and 0.6 for θ as starting values
[param_est_11_ex10, ~,~, cov_11_ex10] = CML(@lik_func_min, @lik_contrib, MA1_realiz_ex10, [10, 0.6, 1]', 2, 1, options);

fprintf('\n\nConditional ML with MA(1):\n');
fprintf('                      μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f \n',[["True          ";"Estimate      "; "Standard Error"], [10; param_est_11_ex10(1); sqrt(cov_11_ex10(1,1))], [5; param_est_11_ex10(2); sqrt(cov_11_ex10(2,2))], [1; param_est_11_ex10(3); sqrt(cov_11_ex10(3,3))]]')

% Interpretation
fprintf('\n- estimate for μ good, but not for θ and σ^2\n');
fprintf('  -> because non-invertible process (θ > 1) !\n');
fprintf('  -> invert it using the twin process with: θ_tilde = 1/θ and σ^2_tilde = θ^2 * σ^2\n');
fprintf('   -> μ already good estimate because same value for both alternatives (μ_tilde = μ)  \n');

% derive parameters of invertible "twin" of MA(1)
theta_tilde = 1/5;
sigma_2_tilde = 5^2 * 1^2;

% Interpretation
fprintf('\n\nCompare conditional ML estimates with parameters of invertible "good twin" MA(1):\n\n');
fprintf('                      μ          θ         σ^2\n');
fprintf('%s %10.3f %10.3f %10.3f\n',[["True (Twin)   ";"Estimate      "; "Standard Error"], [10; param_est_11_ex10(1); sqrt(cov_11_ex10(1,1))], [theta_tilde; param_est_11_ex10(2); sqrt(cov_11_ex10(2,2))], [sigma_2_tilde; param_est_11_ex10(3); sqrt(cov_11_ex10(3,3))]]')

fprintf('\n-> Estimates for θ and σ^2 much better (for μ stays the same) \n');