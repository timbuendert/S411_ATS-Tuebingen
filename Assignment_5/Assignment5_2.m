clear
clc
close all

rng(602); % set seed to ensure replicability

fprintf("Exercise 2 \n")
fprintf("------------------------------- \n")

%% Exercise 2-1

AR1_realiz = func_AR1(1, 1, 0, 0, 1000); % simulate RW with drift
AR1_est = est_OLS(AR1_realiz, 2); % estimation (as case 3 has same estimated model as case 2)
fprintf("\nCase 3:\nρ_hat = %.2f \ns.e.(ρ_hat) = %.2f \nt statistic = %.2f \n", [AR1_est(1), AR1_est(2), AR1_est(4)])


%% Exercise 2-2

n = 10000;
AR_10k_3 = zeros(n, 4);

% generate n realizations of RW with drift and estimate ρ for case 3
for i = 1:n
    AR1_realiz = func_AR1(1, 1, 0, 0, 1000); % simulate RW with drift   
    AR_10k_3(i,:) = est_OLS(AR1_realiz, 2); % estimation (as case 3 has same estimated model as case 2)
end

%% Exercise 2-3

[f3,xi3] = ksdensity(AR_10k_3(:,4)); % estimate kernel density of t statistic
snd = normpdf(-5:0.1:4.9,0,1); % evaluate density of SND based on sequence

figure;
plot(xi3, f3, "b", -5:0.1:4.9, snd, "r"); % plot both densities
title("Comparision of t statistic under H_0 with SND for case 3 ");
legend('t statistic','SND', 'Location', 'northeast')
xlabel('value', 'Fontsize', 14);
ylabel('p', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation:
% - distribution of t statistic almost identical to SND

%% Exercise 2-4
fprintf("\nTrend Stationary process:")

% func_AR1(c, phi, yo, trend, T)
AR1TS_realiz = func_AR1(1, 0.85, (1/(1-0.85)), 0.1, 100); % simulate trend-stationary AR(1)
AR1TS_est = est_OLS(AR1TS_realiz, 2); % estimation (as case 3 has same estimated model as case 2)
fprintf("\nCase 3:\nρ_hat = %.2f \nt statistic = %.2f \n", [AR1TS_est(1), AR1TS_est(4)])

% repeat simulation and computation of test statistic n = 10000 times
n = 10000;
ARTS_10k = zeros(n, 4);

% generate n realizations of trend-stationary process and estimate ρ for case 3
for i = 1:n
    AR1_realiz = func_AR1(1, 0.85, (1/(1-0.85)), 0.1, 100); % simulate trend-stationary process  
    ARTS_10k(i,:) = est_OLS(AR1_realiz, 2); % estimation (as case 3 has same estimated model as case 2)
end

% estimate kernel density of the test statistic under the trend-stationary alternative
[f1,xi1] = ksdensity(ARTS_10k(:,4)); 

% plot kernel densities of t statistics under H_0 and trend-stationary alternative
figure;
plot(xi1, f1, "b", xi3, f3, "r");
title("Comparision of t statistic for case 3 ");
legend('t statistic trend-stationary','t statistic unit root', 'Location', 'northeast')
xlabel('value', 'Fontsize', 14);
ylabel('p', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation:
% - distribution of t statistics almost identical: both centered around 0
% -> just distribution under unit root less peaked

%% Exercise 2-5
% calculate relative frequency with which a Type II error is made (for t
% statistic only) 
[~, t_error] = eval_type2error(AR_10k_3(:,3), AR_10k_3(:,4), ARTS_10k(:,3), ARTS_10k(:,4));
fprintf("\n\nType II error of t statistic = %.2f \nHence, power of the test = %.2f\n", [t_error, 1-t_error]) % only consider t statistic here

%% Exercise 2-6

fprintf("\n----------------------------------------\n")

% set parameters
rho = 0.3; % 0.3, 0.85
T = 1000; % 100, 1000

fprintf("\nρ = %.2f, T = %d:\n", [rho, T])

% repeat simulation and computation of test statistic n = 10000 times
n = 10000;
ARTS_10k = zeros(n, 4);

% generate n realizations of trend-stationary process and estimate ρ for case 3
for i = 1:n
    AR1_realiz = func_AR1(1, rho, (1/(1-rho)), 0.1, T); % simulate trend-stationary process  
    ARTS_10k(i,:) = est_OLS(AR1_realiz, 2); % estimation (as case 3 has same estimated model as case 2)
end

% compute kernel density of t statistic
[f1,xi1] = ksdensity(ARTS_10k(:,4)); 

% plot kernel densities of t statistics under H_0 and trend-stationary alternative
figure;
plot(xi1, f1, "b", xi3, f3, "r");
title("Comparision of t statistics for case 3: ρ = ", rho);
legend('t statistic trend-stationary','t statistic unit root', 'Location', 'northeast')
xlabel('value', 'Fontsize', 14);
ylabel('p', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% calculate relative frequency with which a Type II error is made (for t
% statistic only)
[test_error, t_error] = eval_type2error(AR_10k_3(:,3), AR_10k_3(:,4), ARTS_10k(:,3), ARTS_10k(:,4));
fprintf("\nType II error of t statistic = %.2f \nHence, power of the test = %.2f\n", [t_error, 1-t_error]) % only consider t statistic here


% Interpretation:
% ρ = 0.3 & T = 100: less overlapping as trend stationary distribution 
% shifted to the left (power = 0.09 -> for test statistic: power = 1)
% ρ = 0.3 & T = 1000: trend-stationary disitribution shifted to left and
% a lot more peaked (power = 0)

% ρ = 0.85 & T = 100: see Ex. 2-4 -> distribution of t statistics almost 
% identical (power = 0)
% ρ = 0.85 & T = 1000: as for T = 100 -> almost identical distributions,
% just the trend-stationary one a lot more peaked (power = 0)

% ------------------------------------

% Using estimation with Case 4:
% ρ = 0.3 & T = 100: power = 1
% ρ = 0.3 & T = 1000: power = 1

% ρ = 0.85 & T = 100: power = 1
% ρ = 0.85 & T = 1000: power = 1
