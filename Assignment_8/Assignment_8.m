clear
clc
close all

fprintf("Exercise 1 \n")
fprintf("------------------------------- \n\n")


%% Exercise 1-1

data = xlsread('data_abx.xlsx');
home_CAD = data(:,1);
foreign_USD = data(:,2);
rate_CADUSD = data(:,3);
T = size(data,1);

% first column: time series of the ln home market (TSX) prices
% second column: ln foreign market (NYSE) prices
% third column: ln CAD/USD exchange rate
% fourth column: dummy variable indicating the first observation of a day (indicator = 1: first observation of the day, otherwise zero)

%% Exercise 1-2

foreign_CAD = log(exp(foreign_USD) .* exp(rate_CADUSD));

%% Exercise 1-3

figure
plot(1:T, home_CAD, 1:T, foreign_CAD)
title("Comparision of series of home & foreign market")
legend('home market (TSX)', 'foreign market (NYSE)', 'Location', 'southeast')
xlabel('t', 'Fontsize', 14);
ylabel('CAD', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation:
% - strong co-movements, barely without overlapping
% - no deterministic deviation -> cointegration

%% Exercise 1-4

[~, ~, ~, home_tstat] = DF_test(home_CAD, 2);
[~, ~, ~, foreign_tstat] = DF_test(foreign_CAD, 2);

fprintf("Dickey-Fuller t-statistic for home market: %.2f\n", home_tstat)
fprintf("Dickey-Fuller t-statistic for foreign market: %.2f\n", foreign_tstat)

DF2_crit = -2.86;
fprintf("\nCritical value for DF case 2 (α = 5%%): %.2f\n", DF2_crit)

% For both time series, the HO of a unit root cannot be rejected (at α =
% 5%).

%% Exercise 1-5

% would suggest a = (1 -1)' (from theory of PPP) -> as series very close
% also graph: difference almost 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 2-1

fprintf("\nExercise 2 \n")
fprintf("------------------------------- \n\n")

X = [ones(T,1), foreign_CAD];
y = home_CAD;
b_est = inv(X'*X)*(X'*y);
b1 = b_est(2);
res = y - b_est(1)*ones(T,1) - b1*foreign_CAD;

%% Exercise 2-2
a = [1 -b1];
fprintf("Estimated cointegration vector a = [%.4f %.4f]\n\n", [a(1), a(2)])

% Normalized: first element of a = 1

%% Exercise 2-3
[~, ~, ~, res_tstat] = DF_test(res, 1);

fprintf("Dickey-Fuller t-statistic for residual series: %.2f\n", res_tstat)

DF1_crit = -2.76;
fprintf("Critical value for DF case 1 (α = 5%%): %.2f\n", DF1_crit)

% the HO of a unit root can be rejected (at α = 5%)
% -> so, have case of cointegration: two series for which H0 unit root
% cannot be rejected (non-stationary) + linear combination with a = [1 1]:
% z_t = 1*y1t - b0 - 1*y2t -> reject H0 unit root (stationary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 3-1

fprintf("\nExercise 3 \n")
fprintf("------------------------------- \n\n")

%home_diff = diff(home_CAD);
%foreign_diff = diff(foreign_CAD); see try at very end

%% Exercise 3-2
index = find(data(:,4) == 0);

home_CAD_new = home_CAD(index);
foreign_CAD_new = foreign_CAD(index);
res_new = res(index);

%% Exercise 3-3
home_diff = diff(home_CAD_new);
foreign_diff = diff(foreign_CAD_new);


T_new = size(home_diff, 1) - 1;

const = ones(T_new,1);
home_lag = lagmatrix(home_diff, 1);
foreign_lag = lagmatrix(foreign_diff, 1);

X_new = [const(2:T_new), res_new(3:T_new+1), home_lag(2:T_new), foreign_lag(2:T_new)];
y_new = [home_diff(3:T_new+1), foreign_diff(3:T_new+1)];

b_est_new = inv(X_new'*X_new)*(X_new'*y_new);


%% Exercise 3-4

a10 = b_est_new(1,1);
a20 = b_est_new(1,2);
gamma1 = b_est_new(2,1);
gamma2 = b_est_new(2,2);
a11 = b_est_new(3,1);
a21 = b_est_new(3,2);
a12 = b_est_new(4,1);
a22 = b_est_new(4,2);

fprintf("Δy_1t = %.2f + %.2f * e_t-1 + %.2f * Δy_1t-1 + * %.2f * Δy_2t-1 + u_1t\n", [a10, gamma1, a11, a12])
fprintf("gamma_1 = %.2f\n" , gamma1)

fprintf("\nΔy_2t = %.2f + %.2f * e_t-1 + %.2f * Δy_1t-1 + * %.2f * Δy_2t-1 + u_2t\n", [a20, gamma2, a21, a22])
fprintf("gamma_2 = %.2f\n" , gamma2)


% Interpretation:
% - gamma1 and gamma2 have different signs because they adjust for the
% deviation from the equilibrium (PPP) -> different signs shows that
% mechanism (buying and selling) brings the two prices together
% error correction term


%% other try: first take differences, then exclude specific differences
% -> not yield alternative signs for gammas

% home_diff = [0;diff(home_CAD)];
% foreign_diff = [0;diff(foreign_CAD)];
% 
% df = [home_diff, foreign_diff, res, data(:,4)];
% df_new = df(df(:,4) == 0,:);
% 
% T_new = size(df_new, 1);
% 
% const = ones(T_new,1);
% home_lag = lagmatrix(df_new(:,1), 1);
% foreign_lag = lagmatrix(df_new(:,2), 1);
% res_new = df_new(:,3);
% 
% X_neww = [const(4:T_new-3), res_new(4:T_new-3), home_lag(4:T_new-3), foreign_lag(4:T_new-3)];
% y_neww = [home_diff(5:T_new), foreign_diff(5:T_new)];
% 
% b_est_new = inv(X_new'*X_new)*(X_new'*y_new)
