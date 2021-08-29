clear
clc
close all

rng(602); % set seed to ensure replicability

%% Determine case

case_test = 4; 

%%%%%%%%%%%%%%%%%%%%%%%

if case_test == 1
    fprintf('CASE 1: y_t = ρ * y_t−1 + u_t\n');
elseif case_test == 2
    fprintf('CASE 2: y_t = α + ρ*y_t−1 + u_t\n');
elseif case_test == 4
    fprintf('CASE 4: y_t = α + δ*t + ρ*y_t−1 + u_t\n');
end

%% Estimation Part

% expand for-loop for first parts of exercises
for i = 1:1    
%%%%%%%%%%%%%%%%%%%%%% COMMENT LINE
% %% Exercise 1-1
% T = 100; % set number of observations
% RW = func_RW(0, T); % call function based on parameters
% 
% figure;
% plot(1:T+1, RW, 'k'); % plot realization of process
% title('Realization of Random Walk');
% xlabel('i', 'Fontsize', 14);
% ylabel('y', 'Fontsize', 14);
% set(gca, 'Fontsize', 12);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Exercise 2-1
% fprintf('Exercise 2-1 ------------------------------------------------\n');
% 
% X = RW(1:T);
% 
% 
% if case_test == 1
%     X = RW(1:T);
% elseif case_test == 2
%     X = [ones(T,1), RW(1:T)];
% elseif case_test == 4
%     X = [ones(T,1), (1:1:T)', RW(1:T)];
% end
% 
% y = RW(2:T+1);
% 
% b_est = inv(X'*X)*(X'*y);
% 
% if case_test == 1
%     rho_est = b_est(1);
%     fprintf('\n Estimate for ρ: %.3f\n',rho_est)
% 
% elseif case_test == 2
%     alpha_est = b_est(1);
%     rho_est = b_est(2);
%     fprintf('\nEstimate for α: %.3f\nEstimate for ρ: %.3f\n', [alpha_est, rho_est]')
%     
% elseif case_test == 4
%     alpha_est = b_est(1);
%     delta_est = b_est(2);
%     rho_est = b_est(3);
%     fprintf('\nEstimate for α: %.3f\nEstimate for δ: %.3f\nEstimate for ρ: %.3f\n', [alpha_est, delta_est, rho_est]')
%    
% end
% 
% 
% %% Exercise 2-2
% fprintf('\n\nExercise 2-2 ------------------------------------------------\n');
% 
% if case_test == 1
%     s_sq = (1/(T-1))*sum((RW(2:T+1)-rho_est*RW(1:T)).^2);
%     var_rho = s_sq * inv(X'*X);
%     se_rho = sqrt(var_rho);
% 
% elseif case_test == 2
%     s_sq = (1/(T-2))*sum((RW(2:T+1)-alpha_est-rho_est*(RW(1:T))).^2);
%     var_b = s_sq * inv(X'*X);
%     se_rho = sqrt(var_b(2,2));
%     
% elseif case_test == 4
%     s_sq = (1/(T-3))*sum((RW(2:T+1)-alpha_est-delta_est*(1:1:T)'-rho_est*(RW(1:T))).^2);
%     var_b = s_sq * inv(X'*X);
%     se_rho = sqrt(var_b(3,3));
%    
% end
%     
% t_stat_rho = (rho_est-1)/se_rho;
% 
% fprintf('\n t-statistic: %.3f\n',t_stat_rho)
% 
% 
% %% Exercise 2-3
% fprintf('\n\nExercise 2-3 ------------------------------------------------\n');
% 
% fprintf('\n Test statistic: %.3f\n',T*(rho_est-1))
% 
% 
% % -> see est_OLS.m
%%%%%%%%%%%%%%%%%%%%%% COMMENT LINE
end 
%

n = 10000; % set number of iterations
quant = [0.01 0.025 0.05 0.1 0.9 0.95 0.975 0.99]; % define quantiles

for T = [100, 1000] % vary length of realization (T)

    RW_10k = zeros(n,4); % define empty matrix for saving the estimation results

    for i = 1:n % iterate through n
        RW_real = func_RW(0, T); % simulate RW
        RW_10k(i,:) = est_OLS(RW_real, case_test); % save results from estimation in matrix
    end

    t_quant = quantile(RW_10k(:,4),quant); % compute quantiles of t-statistic
    test_quant = quantile(RW_10k(:,3),quant); % compute quantiles of test statistic

    fprintf('\n\nQuantiles of t- and test statistic for T = %i:\n\n', T);
    fprintf('%20.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',quant')
    fprintf('%s %4.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',["t-statistic   ", t_quant]')
    fprintf('%s %4.2f %7.2f %6.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',["Test-statistic", test_quant]')
end

% checking with Tables B.5 (T (ρˆ − 1)) and B.6 (t−statistic) in Hamilton (1994, p.762+763):
%  -> close enough, correct implementation


% using T = 1000 and n = 10000: visiualize the distributions of the n
% values

% ρ_est − 1
series_1 = RW_10k(:,1) - 1;
[f1,xi1] = ksdensity(series_1); 

figure;
plot(xi1, f1, 'k'); % plot realization of process
title(['Case ', num2str(case_test),': Distribution of ρ(rho est) − 1']);
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% sqrt(T) * ρ_est − 1
series_2 = sqrt(T) * (RW_10k(:,1) - 1);
[f2,xi2] = ksdensity(series_2); 

figure;
plot(xi2, f2, 'k'); % plot realization of process
title(['Case ', num2str(case_test),': Distribution of sqrt(T) * (ρ(rho est) − 1)']);
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% Test-statistic = T * ρ_est − 1
series_3 = RW_10k(:,3);
[f3,xi3] = ksdensity(series_3); 

snd = normpdf(xi3,0,1);

figure;
plot(xi3, f3, xi3, snd, 'k'); % plot realization of process
title(['Case ', num2str(case_test),': Distribution of T * (ρ(rho est) − 1)']);
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
legend('Test-statistic','N(0,1)')
set(gca, 'Fontsize', 12);

% t-statistic
[f4,xi4] = ksdensity(RW_10k(:,4)); 

snd = normpdf(xi4,0,1);

figure;
plot(xi4, f4, xi4, snd, 'k'); % plot realization of process
title(['Case ', num2str(case_test),': Distribution of t-statistic']);
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
legend('t-statistic','N(0,1)')
set(gca, 'Fontsize', 12);


fprintf('\nSee graphs\n');


%% Interpretation

%%% Case 1:
% ρ_est − 1: very centered around 0, almost degenerate
% sqrt(T) * ρ_est − 1: less, but still very centered around 0
% Test-statistic = T * ρ_est − 1: still quite centered around 0 but peaked half as much as N(0,1) -> not close to degeneration
% t-statistic: very similar to N(0,1), but slightly more peaked and shifted to the left

%%% Case 2:
% ρ_est − 1: very centered around 0, almost degenerate
% sqrt(T) * ρ_est − 1: less, but still very centered around 0
% Test-statistic = T * ρ_est − 1: way less peaked than N(0,1) and shifted to the left -> not close to degeneration
% t-statistic: very similar to N(0,1), but slightly more peaked and shifted to the left

%%% Case 4:
% ρ_est − 1: very centered around 0, almost degenerate
% sqrt(T) * ρ_est − 1: less, but still very centered around 0
% Test-statistic = T * ρ_est − 1: even less peaked than Case 3 and more shifted to the left -> not close to degeneration
% t-statistic: very similar to N(0,1), but slightly more peaked and shifted to the left
