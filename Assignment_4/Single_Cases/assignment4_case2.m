clear
clc
close all

rng(603); % set seed to ensure replicability

fprintf('CASE 2: y_t = α + ρ*y_t−1 + u_t\n\n');

%% Exercise 3-1
T = 100; % set number of observations
RW = func_RW(0, T); % call function based on parameters

figure;
plot(1:T+1, RW, 'k'); % plot realization of process
title('Realization of Random Walk');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 3-2
fprintf('Exercise 3-2 ------------------------------------------------\n');

X = [ones(T,1), RW(1:T)];
y = RW(2:T+1);

b_est = inv(X'*X)*(X'*y);
alpha_est = b_est(1);
rho_est = b_est(2);
fprintf('\nEstimate for α: %.3f\nEstimate for ρ: %.3f\n', [alpha_est, rho_est]')

%% Exercise 3-3
fprintf('\n\nExercise 3-3 ------------------------------------------------\n');

s_sq = (1/(T-2))*sum((RW(2:T+1)-alpha_est-rho_est*(RW(1:T))).^2);
var_b = s_sq * inv(X'*X);
se_rho = sqrt(var_b(2,2));
t_stat_rho = (rho_est-1)/se_rho;

fprintf('\n t-statistic: %.3f\n',t_stat_rho)

%% Exercise 3-4
fprintf('\n\nExercise 3-4 ------------------------------------------------\n');


fprintf('\n Test statistic: %.3f\n',T*(rho_est-1))


% -> see est_OLS2.m


n = 10000;
quant = [0.01 0.025 0.05 0.1 0.9 0.95 0.975 0.99];

for T = [100, 1000]

    RW_10k = zeros(n,4);

    for i = 1:n
        RW_real = func_RW(0, T);
        RW_10k(i,:) = est_OLS2(RW_real);
    end

    t_quant = quantile(RW_10k(:,4),quant);
    test_quant = quantile(RW_10k(:,3),quant);

    fprintf('\n\nQuantiles of t- and test statistic for T = %i:\n\n', T);
    fprintf('%20.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',quant')
    fprintf('%s %4.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',["t-statistic   ", t_quant]')
    fprintf('%s %4.2f %7.2f %6.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',["Test-statistic", test_quant]')
end

% checking with Tables B.5 (T (ρˆ − 1)) and B.6 (t−statistic) in Hamilton (1994, p.762+763):
%  -> close enough, correct implementation


% using T = 1000 and n = 10000

% ρ_est − 1
series_1 = RW_10k(:,1) - 1;
[f1,xi1] = ksdensity(series_1); 

figure;
plot(xi1, f1, 'k'); % plot realization of process
title('Case 2: Distribution of ρ(est) − 1');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% sqrt(T) * ρ_est − 1
series_2 = sqrt(T) * (RW_10k(:,1) - 1);
[f2,xi2] = ksdensity(series_2); 

figure;
plot(xi2, f2, 'k'); % plot realization of process
title('Case 2: Distribution of sqrt(T) * ρ(est) − 1');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% Test-statistic = T * ρ_est − 1
series_3 = RW_10k(:,3);
[f3,xi3] = ksdensity(series_3); 

figure;
plot(xi3, f3, 'k'); % plot realization of process
title('Case 2: Distribution of T * ρ(est) − 1');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% t-statistic
[f4,xi4] = ksdensity(RW_10k(:,4)); 

figure;
plot(xi4, f4, 'k'); % plot realization of process
title('Case 2: Distribution of t-statistic');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


fprintf('\nSee graphs\n');

% Interpretation:
% -
