clear
clc
close all

rng(101); % set seed to ensure replicability

T = 100; % set number of observations


%% Exercise 1
% -> see func_MA2()


%% Exercise 2
MAa = func_MA2(0,0.7,0.3,T); % call function based on parameters in a)
MAb = func_MA2(4,5,-17,T); % call function based on parameters in b)


%% Exercise 3
x = (1:T)'; % create corresponding sequence for plotting

% a)
figure;
plot(x, MAa, x, ones(T,1)*0, 'k'); % plot realization and expected value
title('Realization of MA(2) in part a)');
legend('MA(2)', 'E(Y)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% b)
figure;
plot(x,MAb,  x, ones(T,1)*4, 'k'); % plot realization and expected value
title('Realization of MA(2) in part b)');
legend('MA(2)', 'E(Y)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);


%% Exercise 4
% a)
MAa_ac = (autocorr(MAa)); % compute empirical autocorrelations
disp(['The first-order empirical autocorrelation of the MA(2) of a) is: ' num2str(MAa_ac(2))]); % select autocorr. of first order
disp(['The first-order theoretical autocorrelation of the MA(2) of a) is: ' num2str((1*(0.7+0.3*0.7)) / (1+0.7^2+0.3^2))]);
disp(['The deviation of the MA(2) of a) is: ' num2str(100*(1-(MAa_ac(2) / ((1*(0.7+0.3*0.7)) / (1+0.7^2+0.3^2))))) '%']);

% b) -> see comments for 4.a)
MAb_ac = (autocorr(MAb));
disp(['The first-order empirical autocorrelation of the MA(2) of b) is: ' num2str(MAb_ac(2))]);
disp(['The first-order theoretical autocorrelation of the MA(2) of a) is: ' num2str((1*(5+5*(-17)))/ (1+5^2+(-17)^2))]);
disp(['The deviation of the MA(2) of a) is: ' num2str(100*(1-(MAb_ac(2) / ((1*(5+5*(-17))) / (1+5^2+(-17)^2))))) '%']);

% Interpretation
% - not too much of an absolute difference between empirical and theoretical autocorrelations 
% -> close considering T = 100


%% Exercise 5
n_ensembles = 30;

% a)
MAa_m = zeros(T, n_ensembles); % initialize matrix for ensembles
for i = 1:n_ensembles
    MAa_m(:,i) = func_MA2(0,0.7,0.3,T); % iteratively add realization with parameters of a)
end 
MAa_ensMeans = mean(MAa_m, 2); % compute ensemble means across T observations
MAa_ensVar = var(MAa_m, 0, 2); % compute ensemble variances across T observations

figure;
plot(x, MAa_ensMeans,x, MAa_ensVar, 'r', x, ones(T,1)*0, 'k', x, ones(T,1)*((0.7)^2 + (0.3)^2 + 1), 'k'); % plot ensemble means and variances + theoretical values
title('Ensemble means and variances of MA(2) in part a)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% b) -> see comments for 5.a)
MAb_m = zeros(T, n_ensembles);
for i = 1:n_ensembles
    MAb_m(:,i) = func_MA2(4,5,-17,T);
end 
MAb_ensMeans = mean(MAb_m, 2);
MAb_ensVar = var(MAb_m, 0, 2);

figure;
plot(x, MAb_ensMeans, x, MAb_ensVar, 'r', x, ones(T,1)*4, 'k', x, ones(T,1)*(5^2 + (-17)^2 + 1), 'k');
title('Ensemble means and variances of MA(2) in part b)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% as n_ensembles -> infinity: means converge to theoretical values

%% Exercise 6
% a)
MAa_corr = corr(MAa_m'); % compute autocorrelation matrix

% initialize autocorr. matrix -> T-1 because there is no first-order correlation between y1 and y100 -> one missing
MAa_corr_ext = zeros(T-1,1);
for i = 1:(size(MAa_m,1)-1) % iterate through autocorrelations
    MAa_corr_ext(i,1) = MAa_corr(i, i+1); % read out selected first-order autocorrelations
end 

figure;
plot(x(1:T-1), MAa_corr_ext, x(1:T-1), ones(T-1,1)*((0.7+0.3*0.7) / (1+0.7^2+0.3^2)), 'k'); % plot first-order autocorrelations + theoretical values
title('Autocorrelations of ensembles of MA(2) in part a)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - Autocorrelations vary greatly over the 100 observations
% -> mean: roughly 0.55 (values vary around theoretical value)
% - a): only positive -> positive lines relationship between observations with one period in-between


% b) -> see commments for 6.a)
MAb_corr = corr(MAb_m');
MAb_corr_ext = zeros(T-1,1);
for i = 1:(size(MAa_m,1)-1)
    MAb_corr_ext(i,1) = MAb_corr(i, i+1);
end 

figure;
plot(x(1:T-1), MAb_corr_ext, x(1:T-1), ones(T-1,1)*((5+5*(-17))/ (1+5^2+(-17)^2)), 'k');
title('Autocorrelations of ensembles of MA(2) in part b)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - Autocorrelations vary greatly over the 100 observations
% -> mean: roughly -0.3 (values vary around theoretical value)
% - only negative -> negative linear relationship between observations with one period in-between


%% Exercise 7
% Stationarity does not depend on phi1 and ph2: all MA(q) process are stationary

% Process a)
% - process looks like distribution does not change over time (1.3)
% - ensemble means (1.5) fluctuate a lot but around a constant value (roughly 1.5)
% - first-order autocorrelations fluctuate around 0.5 (1.6) -> < 1
% -> stationary (empirical values fluctuate around theoretical value -> as ensembles and T grow: more towards theory)

% Process b)
% - process looks like distribution does not change over time: very erratic, higher absolute values (1.3)
% - ensemble means (1.5) fluctuate a lot but around a constant value (roughly 300) -> higher fluctuations
% - first-order autocorrelations fluctuate around -0.3 (1.6) -> < 1
% -> stationary (empirical values fluctuate around theoretical value -> as ensembles and T grow: more towards theory)