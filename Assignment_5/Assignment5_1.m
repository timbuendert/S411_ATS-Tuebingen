clear
clc
close all

fprintf("Exercise 1 \n")
fprintf("------------------------------- \n")

%% Exercise 1-1
rng(602); % set seed to ensure replicability

fprintf('CASE 1: y_t = ρ * y_t−1 + u_t\n');
fprintf('CASE 2: y_t = α + ρ*y_t−1 + u_t\n');
fprintf('CASE 4: y_t = α + δ*t + ρ*y_t−1 + u_t\n');

T = 1000;
n = 10000;

RW_10k_1 = zeros(n, 4);
RW_10k_2 = zeros(n, 4);
RW_10k_4 = zeros(n, 4);

% generate n realizations of RW and estimate ρ for cases 1, 2, 4

for i = 1:n
    RW_real = func_RW(0, T); % simulate RW     
    RW_10k_1(i,:) = est_OLS(RW_real, 1); % case 1
    RW_10k_2(i,:) = est_OLS(RW_real, 2); % case 2
    RW_10k_4(i,:) = est_OLS(RW_real, 4); % case 4
end

% generate kernel densities and mean of ρ_est
[f1,xi1] = ksdensity(RW_10k_1(:,1)); 
mean_1 = mean(RW_10k_1(:,1));

[f2,xi2] = ksdensity(RW_10k_2(:,1)); 
mean_2 = mean(RW_10k_2(:,1));

[f4,xi4] = ksdensity(RW_10k_4(:,1)); 
mean_4 = mean(RW_10k_4(:,1));

% plot of the three kernel densities + mean
figure;
plot(xi1, f1, "b", xi2, f2, "g", xi4, f4, "r"); % plot realization of process
xline(mean_1, "--b");
xline(mean_2, "--g");
xline(mean_4, "--r");
title("Comparision of ρ_{est} for cases 1, 2, 4 ");
legend('Case 1','Case 2','Case 4', 'Location', 'northwest')
xlabel('ρ', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation:
% - all center around 1, while case 1 almost degenerate 
%   -> case 2 less and case 4 even less (towards 0.99 and less peaked)
% -> case 1 most accurately estimates ρ (reasonable as case 2 expects
%    intercept and case 4 expects intercept + time trend)

%% Exercise 1-2
rng(602); % set seed to ensure replicability

% plot distributions of test statistics for cases 1, 2, 4 in one graph:

% Test statistic
[test1,testx1] =  ksdensity(RW_10k_1(:,3));  % Case 1
[test2,testx2] =  ksdensity(RW_10k_2(:,3));  % Case 2
[test4,testx4] =  ksdensity(RW_10k_4(:,3));  % Case 4

figure;
plot(testx1, test1, "b", testx2, test2, "g", testx4, test4, "r"); % plot realization of process
title("Comparision of test-statistic for cases 1, 2, 4 ");
legend('Case 1','Case 2','Case 4', 'Location', 'northwest')
xlabel('test-statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% t statistic
[t1,tx1] = ksdensity(RW_10k_1(:,4));  % Case 1
[t2,tx2] = ksdensity(RW_10k_2(:,4));  % Case 2
[t4,tx4] = ksdensity(RW_10k_4(:,4));  % Case 4

figure;
plot(tx1, t1, "b", tx2, t2, "g", tx4, t4, "r"); % plot realization of process
title("Comparision of t-statistic for cases 1, 2, 4 ");
legend('Case 1','Case 2','Case 4', 'Location', 'northwest')
xlabel('t-statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation:
% - both statistics of case 1 center around 0 while case 2 left slightly 
% to the left and case 4 even more to the left
% -> as before: reasonable given that case 2 expects an intercept and
% case 4 an intercept + time trend which are not present in the
% data - generating process

%% Exercise 1-3
rng(602); % set seed to ensure replicability

fprintf("\n----------------------------------------\n")
fprintf("\nρ = 0.85:\n")

% Sample calculation: 
%AR1_realiz_st = func_AR1(1, 0.85, (1/(1-0.85)), 0, 100);
%ouput_st = est_OLS(AR1_realiz_st, 1);
%fprintf("\nSample calculation for stationary AR(1):\n")
%fprintf("Test statistic: %.2f // T-statistic: %.2f \n\n", [ouput_st(3), ouput_st(4)])

T = 100;
n = 10000;

RW_10k_1_st = zeros(n, 4);
RW_10k_2_st = zeros(n, 4);
RW_10k_4_st = zeros(n, 4);

% generate n realizations of stationary AR(1) and estimate parameters
% -> compute test statistics

for i = 1:n
    RW_real = func_AR1(1, 0.85, (1/1-0.85), 0, T); % simulate stationary AR(1): starting value important!     
    RW_10k_1_st(i,:) = est_OLS(RW_real, 1); % case 1
    RW_10k_2_st(i,:) = est_OLS(RW_real, 2); % case 2
    RW_10k_4_st(i,:) = est_OLS(RW_real, 4); % case 4
end

% plot stationary and RW distributions of test statistics in one graph:

% Case 1
[test1_st,testx1_st] =  ksdensity(RW_10k_1_st(:,3));  % distribution of test statistic
[t1_st,tx1_st] = ksdensity(RW_10k_1_st(:,4));  % distribution of t statistic

figure;
plot(testx1, test1, "b", tx1, t1, "g", testx1_st, test1_st, "r", tx1_st, t1_st, "y");
title("Comparision of test statistics for case 1 (ρ = 0.85)");
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% Case 2
[test2_st,testx2_st] =  ksdensity(RW_10k_2_st(:,3));  % distribution of test statistic
[t2_st,tx2_st] = ksdensity(RW_10k_2_st(:,4));  % distribution of t statistic

figure;
plot(testx2, test2, "b", tx2, t2, "g", testx2_st, test2_st, "r", tx2_st, t2_st, "y");
title("Comparision of test statistics for case 2 (ρ = 0.85)");
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Case 4
[test4_st,testx4_st] =  ksdensity(RW_10k_4_st(:,3));  % distribution of test statistic
[t4_st,tx4_st] = ksdensity(RW_10k_4_st(:,4));  % distribution of t statistic

figure;
plot(testx4, test4, "b", tx4, t4, "g", testx4_st, test4_st, "r", tx4_st, t4_st, "y");
title("Comparision of test statistics for case 4 (ρ = 0.85)");
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

%% Exercise 1-4
% see eval_type2error.m
% type II error probability: probability of failing to reject the incorrect 
% null hypothesis of a unit root

%% Exercise 1-5

% evaluate type II errors based on RW and stationary statistics

[test_error1, t_error1] = eval_type2error(RW_10k_1(:,3), RW_10k_1(:,4), RW_10k_1_st(:,3), RW_10k_1_st(:,4)); % case 1
fprintf("\nCase 1:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error1, t_error1])

[test_error2, t_error2] = eval_type2error(RW_10k_2(:,3), RW_10k_2(:,4), RW_10k_2_st(:,3), RW_10k_2_st(:,4)); % case 2
fprintf("\nCase 2:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error2, t_error2])

[test_error4, t_error4] = eval_type2error(RW_10k_4(:,3), RW_10k_4(:,4), RW_10k_4_st(:,3), RW_10k_4_st(:,4)); % case 4
fprintf("\nCase 4:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error4, t_error4])

% Interpretation:
% (relative frequency of making a type II error = 0 -> power = 1)

% - Case 1: distributions of statistics for RW and stationary series 
%   overlap greatly + Type II error for t- and test statistic = 1

% - Case 2: less overlapping of disitributions and significantly smaller
%   Type II errors

% - Case 4: overlapping of disitributions and test statistics between case 
%   1 and 2 and significantly smaller

% -> for this process: test 1 < test 4 < test 2 (largest power)
%   -> reasonable given true stationary AR(1) has intercept, so has
%   estimated model of case 2
% -> also, test statistic < t statistic
%   -> best: case 2 + t statistic


%% Exercise 1-6
rng(602); % set seed to ensure replicability

fprintf("\n----------------------------------------\n")

% set parameters
rho = 0.3; % 0.3, 0.85, 0.98
T = 1000; % 100, 1000

fprintf("\nρ = %.2f, T = %d:\n", [rho, T])

n = 10000;

RW_10k_1_st = zeros(n, 4);
RW_10k_2_st = zeros(n, 4);
RW_10k_4_st = zeros(n, 4);

% generate n realizations of stationary AR(1) and estimate parameters
% -> compute test statistics

for i = 1:n
    RW_real = func_AR1(1, rho, (1/(1-rho)), 0, T); % simulate stationary AR(1)     
    RW_10k_1_st(i,:) = est_OLS(RW_real, 1); % case 1
    RW_10k_2_st(i,:) = est_OLS(RW_real, 2); % case 2
    RW_10k_4_st(i,:) = est_OLS(RW_real, 4); % case 4
end


% plot stationary and RW distributions of test statistics in one graph:

% Case 1
[test1_st,testx1_st] =  ksdensity(RW_10k_1_st(:,3));  % distribution of test statistic
[t1_st,tx1_st] = ksdensity(RW_10k_1_st(:,4));  % distribution of t statistic

figure;
plot(testx1, test1, "b", tx1, t1, "g", testx1_st, test1_st, "r", tx1_st, t1_st, "y"); 
title("Comparision of test statistics for case 1 with ρ = ", rho)
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% Case 2
[test2_st,testx2_st] =  ksdensity(RW_10k_2_st(:,3));  % distribution of test statistic
[t2_st,tx2_st] = ksdensity(RW_10k_2_st(:,4));  % distribution of t statistic

figure;
plot(testx2, test2, "b", tx2, t2, "g", testx2_st, test2_st, "r", tx2_st, t2_st, "y");
title("Comparision of test statistics for case 2 with ρ = ", rho)
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Case 4
[test4_st,testx4_st] =  ksdensity(RW_10k_4_st(:,3));  % distribution of test statistic
[t4_st,tx4_st] = ksdensity(RW_10k_4_st(:,4));  % distribution of t statistic

figure;
plot(testx4, test4, "b", tx4, t4, "g", testx4_st, test4_st, "r", tx4_st, t4_st, "y");
title("Comparision of test statistics for case 4 with ρ = ", rho)
legend('Test statistic Unit Root', 't statistic Unit root', 'Test statistic Stationary','t statistic Stationary', 'Location', 'northwest')
xlabel('statistic', 'Fontsize', 14);
ylabel('n', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% evaluate type II errors based on RW and stationary statistics

[test_error1, t_error1] = eval_type2error(RW_10k_1(:,3), RW_10k_1(:,4), RW_10k_1_st(:,3), RW_10k_1_st(:,4)); % case 1
fprintf("\nCase 1:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error1, t_error1])

[test_error2, t_error2] = eval_type2error(RW_10k_2(:,3), RW_10k_2(:,4), RW_10k_2_st(:,3), RW_10k_2_st(:,4)); % case 2
fprintf("\nCase 2:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error2, t_error2])

[test_error4, t_error4] = eval_type2error(RW_10k_4(:,3), RW_10k_4(:,4), RW_10k_4_st(:,3), RW_10k_4_st(:,4)); % case 4
fprintf("\nCase 4:\nType II error of test statistic: %.2f \nType II error of t statistic: %.2f \n", [test_error4, t_error4])

% Interpretation:

% ρ = 0.3 & T = 100: all tests & cases can be chosen (all have power = 1)
% -> case 4: distributions farthest away (case 2 second)   
% ρ = 0.3 & T = 1000: all tests & cases can be chosen (all have power = 1)
% -> case 4: distributions farthest away (case 2 second)

% ρ = 0.85 & T = 100: choose case 2 (see Ex. 1-5) + here, test statistic
% ρ = 0.85 & T = 1000: all tests & cases can be chosen (all have power = 1)
% -> case 4: distributions farthest away (case 2 second)

% ρ = 0.98 & T = 100: choose case 2 + test statistic (but all Type II
% errors > 0.90)
% ρ = 0.98 & T = 1000:  choose case 2 + test statistic (case 4 second and
% case 1 has type II error = 1 for both statistics)
