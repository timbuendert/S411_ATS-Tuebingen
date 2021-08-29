clear
clc
close all

rng(102) % set seed to ensure replicability

T = 100; % set number of observations


%% Exercise 1
% -> see func_AR1()


%% Exercise 2
ARa = func_AR1(0,1,0,T); % call function based on parameters in a) (random walk)
ARb = func_AR1(0.7,1,0,T); % call function based on parameters in b) (random walk with drift)
ARc = func_AR1(4.5,0.1,5,T); % call function based on parameters in c)
ARd = func_AR1(4.4,-0.1,4,T); % call function based on parameters in d)
ARe = func_AR1(4.5,0.9,45,T); % call function based on parameters in e)


%% Exercise 3
x = (1:T)'; % create corresponding sequence for plotting

% a)
figure;
plot(x, ARa, 'k'); % plot realization of process a)
title('Realization of AR(1) in part a)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% b)
figure;
plot(x, ARb, 'k');
title('Realization of AR(1) in part b)'); % plot realization of process b)
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% c)
figure;
plot(x, ARc, x, ones(T,1)*(4.5/(1-0.1)), 'k'); % plot realization of process c) + expected value
title('Realization of AR(1) in part c)');
legend('AR(1)', 'E(Y)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% d)
figure;
plot(x, ARd, x, ones(T,1)*(4.4/(1+0.1)), 'k'); % plot realization of process d) + expected value
title('Realization of AR(1) in part d)');
legend('AR(1)', 'E(Y)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% e)
figure;
plot(x, ARe, x, ones(T,1)*(4.5/(1-0.9)), 'k'); % plot realization of process e) + expected value
title('Realization of AR(1) in part e)');
legend('AR(1)', 'E(Y)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


%% Exercise 4
% - smaller absolute phi: more irratic realization because innovation has more weight; less weight on previous value
% -> also more fluctuations around the mean
% - effect of sign: if < 0 -> oppposite direction -> more fluctuation ; if > 0: same direction

%% Exercise 5
% a)
ARa_ac = (autocorr(ARa)); % create autocorrlation matrix
disp(['The first-order empirical autocorrelation of the AR(1) of a) is: ' num2str(ARa_ac(2))]); % display empirical first-order autocorrelation
disp(['The first-order theoretical autocorrelation of the AR(1) of a) is: ' num2str(1)]); % display theoretical first-order autocorrelation
disp(['The deviation of the AR(1) of a) is: ' num2str(100*(1-(ARa_ac(2)/(1)))) '%']);

% b) -> see comments of 5.a)
ARb_ac = (autocorr(ARb));
disp(['The first-order empirical autocorrelation of the AR(1) of b) is: ' num2str(ARb_ac(2))]);
disp(['The first-order theoretical autocorrelation of the AR(1) of b) is: ' num2str(1)]);
disp(['The deviation of the AR(1) of b) is: ' num2str(100*(1-(ARb_ac(2)/(1)))) '%']);

% c) -> see comments of 5.a)
ARc_ac = (autocorr(ARc));
disp(['The first-order empirical autocorrelation of the AR(1) of c) is: ' num2str(ARc_ac(2))]);
disp(['The first-order theoretical autocorrelation of the AR(1) of c) is: ' num2str(0.1)]);
disp(['The deviation of the AR(1) of c) is: ' num2str(100*(1-(ARc_ac(2)/(0.1)))) '%']);

% d) -> see comments of 5.a)
ARd_ac = (autocorr(ARd));
disp(['The first-order empirical autocorrelation of the AR(1) of d) is: ' num2str(ARd_ac(2))]);
disp(['The first-order theoretical autocorrelation of the AR(1) of d) is: ' num2str(-0.1)]);
disp(['The deviation of the AR(1) of d) is: ' num2str(100*(1-(ARd_ac(2)/(-0.1)))) '%']);

% e) -> see comments of 5.a)
ARe_ac = (autocorr(ARe));
disp(['The first-order empirical autocorrelation of the AR(1) of e) is: ' num2str(ARe_ac(2))]);
disp(['The first-order theoretical autocorrelation of the AR(1) of e) is: ' num2str(0.9)]);
disp(['The deviation of the AR(1) of e) is: ' num2str(100*(1-(ARe_ac(2)/(0.9)))) '%']);

% Interpretation
% - empirical autocorrelations are quite close in absolute terms to theoretical ones considering T = 100


%% Exercise 6
n_ensembles = 30;

% a)
ARa_m = zeros(T, n_ensembles); % initialize ensemble matrix
for i = 1:n_ensembles % iterate through number of ensembles (i.e. 30)
    ARa_m(:,i) = func_AR1(0,1,0,T); % generate new realization with parameters of a)
end 
ARa_ensMeans = mean(ARa_m, 2); % compute ensemble means across T observations
ARa_ensVar = var(ARa_m, 0, 2); % compute ensemble variances across T observations

figure;
plot(x, ARa_ensMeans, x, ARa_ensVar, 'r'); % plot ensemble means and variances
title('Ensemble means and variances of AR(1) in part a)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - similar ro realization: ensemble means grow as t increases
% - ensemble variances remain rather constant
% -> distribution changes over time
% - cannot plot theoretical values as phi = 1


% b) -> see comments of 6.a.)
ARb_m = zeros(T, n_ensembles);
for i = 1:n_ensembles
    ARb_m(:,i) = func_AR1(0.7,1,0,T);
end 
ARb_ensMeans = mean(ARb_m, 2);
ARb_ensVar = var(ARb_m, 0, 2);

figure;
plot(x, ARb_ensMeans, x, ARb_ensVar, 'r'); % plot ensemble means and variances
title('Ensemble means and variances of AR(1) in part b)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - both ensemble means and variances increase with t
% -> distribution changes over time
% - cannot plot theoretical values as phi = 1


% c) -> see comments of 6.a.)
ARc_m = zeros(T, n_ensembles);
for i = 1:n_ensembles
    ARc_m(:,i) = func_AR1(4.5,0.1,5,T);
end 
ARc_ensMeans = mean(ARc_m, 2);
ARc_ensVar = var(ARc_m, 0, 2);

figure;
plot(x, ARc_ensMeans, x, ARc_ensVar, 'r', x, ones(T,1)*(4.5/(1-0.1)), 'k', x, ones(T,1)*(1/(1-(0.1)^2)), 'k'); % plot ensemble means and variances + theoretical values
title('Ensemble means and variances of AR(1) in part c)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - ensemble means and variances fluctuate around theoretical values
% - ensemble means and variances remain roughly constant as t increases
% - as n_ensembles -> infinity: means converge to theoretical values


% d) -> see comments of 6.a.)
ARd_m = zeros(T, n_ensembles);
for i = 1:n_ensembles
    ARd_m(:,i) = func_AR1(4.4,-0.1,4,T);
end 
ARd_ensMeans = mean(ARd_m, 2);
ARd_ensVar = var(ARd_m, 0, 2);

figure;
plot(x, ARd_ensMeans, x, ARd_ensVar, 'r', x, ones(T,1)*(4.4/(1+0.1)), 'k', x, ones(T,1)*(1/(1-(0.1)^2)), 'k'); % plot ensemble means and variances + theoretical values
title('Ensemble means and variances of AR(1) in part d)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - ensemble means and variances fluctuate around theoretical values
% - ensemble means and variances remain roughly constant as t increases
% - as n_ensembles -> infinity: means converge to theoretical values


% e) -> see comments of 6.a.)
ARe_m = zeros(T, n_ensembles);
for i = 1:n_ensembles
    ARe_m(:,i) = func_AR1(4.5,0.9,45,T);
end 
ARe_ensMeans = mean(ARe_m, 2);
ARe_ensVar = var(ARe_m, 0, 2);

figure;
plot(x, ARe_ensMeans, x, ARe_ensVar, 'r', x, ones(T,1)*(4.5/(1-0.9)), 'k', x, ones(T,1)*(1/(1-(0.9)^2)), 'k'); % plot ensemble means and variances + theoretical values
title('Ensemble means and variances of AR(1) in part e)');
legend('Ens. Means', 'Ens. Var.');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% Interpretation
% - ensemble means and variances fluctuate around theoretical values
% - ensemble means and variances remain roughly constant as t increases
% - as n_ensembles -> infinity: means converge to theoretical values
% -> variances converge to theoretical value
% -> because of y0 -> var at t = 0: almost 1 (as innovation small compared to large yo and large phi) -> needs time to forget y0
% (in practice: take T = 150 and cut first 50 t's)

%% Exercise 7
% a)
ARa_corr = corr(ARa_m'); % compute autocorrelation matrix

% initialize autocorr. matrix -> T-1 because there is no first-order correlation between y1 and y100 -> one missing
ARa_corr_ext = zeros(T-1,1);
for i = 1:(size(ARa_m,1)-1) % iterate through autocorrelations
    ARa_corr_ext(i,1) = ARa_corr(i, i+1); % read out selected first-order autocorrelations
end 

figure;
plot(x(1:T-1), ARa_corr_ext, x(1:T-1), ones(T-1,1)*(1), 'k'); % plot first-order autocorrelations + theoretical value
title('Autocorrelations of ensembles of AR(1) in part a)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - first-order autocorrelations converge to 1 (theoretical value) as t increases
% -> strong dependency
% - as n_ensembles -> infinity: autocorrelations converge to theoretical values


% b) -> see comments of 7.a.)
ARb_corr = corr(ARb_m');
ARb_corr_ext = zeros(T-1,1);
for i = 1:(size(ARb_m,1)-1)
    ARb_corr_ext(i,1) = ARb_corr(i, i+1);
end 

figure;
plot(x(1:T-1), ARb_corr_ext, x(1:T-1), ones(T-1,1)*(1), 'k'); % plot first-order autocorrelations + theoretical value
title('Autocorrelations of ensembles of AR(1) in part b)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - first-order autocorrelations converge to 1 (theoretical value) as t increases
% -> strong dependency
% - as n_ensembles -> infinity: autocorrelations converge to theoretical values


% c) -> see comments of 7.a.)
ARc_corr = corr(ARc_m');
ARc_corr_ext = zeros(T-1,1);
for i = 1:(size(ARc_m,1)-1)
    ARc_corr_ext(i,1) = ARc_corr(i, i+1);
end 

figure;
plot(x(1:T-1), ARc_corr_ext, x(1:T-1), ones(T-1,1)*(0.1), 'k'); % plot first-order autocorrelations + theoretical value
title('Autocorrelations of ensembles of AR(1) in part c)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - first-order autocorrelations fluctuate around theoretical value as t increases


% d) -> see comments of 7.a.)
ARd_corr = corr(ARd_m');
ARd_corr_ext = zeros(T-1,1);
for i = 1:(size(ARd_m,1)-1)
    ARd_corr_ext(i,1) = ARd_corr(i, i+1);
end 

figure;
plot(x(1:T-1), ARd_corr_ext, x(1:T-1), ones(T-1,1)*(-0.1), 'k'); % plot first-order autocorrelations + theoretical value
title('Autocorrelations of ensembles of AR(1) in part d)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - first-order autocorrelations fluctuate around theoretical value as t increases


% e) -> see comments of 7.a.)
ARe_corr = corr(ARe_m');
ARe_corr_ext = zeros(T-1,1);
for i = 1:(size(ARe_m,1)-1)
    ARe_corr_ext(i,1) = ARe_corr(i, i+1);
end 

figure;
plot(x(1:T-1), ARe_corr_ext, x(1:T-1), ones(T-1,1)*(0.9), 'k'); % plot first-order autocorrelations + theoretical value
title('Autocorrelations of ensembles of AR(1) in part e)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - first-order autocorrelations converge to theoretical value as t increases (with slightly erratic shape)
% - as n_ensembles -> infinity: autocorrelations converge to theoretical values


%% Exercise 8

% AR(1) stationary: if absolute value of phi < 1
% -> ensemble & time series values converge to theoretical counterparts as t / ensembles -> infinity

% Process a)
% - process looks like distribution does change over time (2.3)
% -> ensemble means (2.6) increase over time 
% - first-order autocorrelations converge to 1 as t increases (2.7)
% -> non-stationary

% Process b)
% - process looks like distribution changes over time (2.3) -> monotonoically increasing values
% -> ensemble means and variances (2.6) increase over time 
% - first-order autocorrelations converge to 1 as t increases (2.7)
% -> non-stationary

% Process c)
% - process looks like distribution does not change over time (2.3)
% - ensemble means (2.6) remain rather constant, same for ensemble variances
% - first-order autocorrelations fluctuate around 0 as t increases (rather positiv) (2.7)
% -> stationary

% Process d)
% - process looks like distribution does not change over time (2.3)
% - ensemble means (2.6) remain rather constant, same for ensemble variances
% - first-order autocorrelations fluctuate around 0 as t increases (rather negative) (2.7)
% -> stationary

% Process e)
% - process looks like distribution does not change too much over time (2.3)
% - ensemble means (2.6) remain rather constant, same for ensemble variances
% - first-order autocorrelations converge to 0.9 as t increases (2.7)
% -> stationary