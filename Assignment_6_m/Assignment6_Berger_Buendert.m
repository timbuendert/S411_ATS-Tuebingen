% Assignment #6 by 
% Leonard Berger (Student-ID: 5629475) 
% Tim-Moritz Bündert (Student-ID: 5635975)

% The following function files accompany this code: 
% - ADF.m
% - adv_ADF.m

% For the code to run smoothly, please add the "VAR2_toolbox" and
% "SWISS_1976_2014.xlsx" to the working directory.

clear
clc
close all


%% Exercise 1-1

% load the data in excel-format into Matlab
SWISS = xlsread("SWISS_1976_2014.xlsx");


%% Exercise 1-2
fprintf('Exercise 1-2 (Section 2 in paper) ------------------------------------------------');


% plot quarterly interest rate (first column of the data)
figure('Position', [0 0 650 450]) % fix size of the plot
plot(SWISS(:,1), "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('\%', 'Fontsize', 14,'interpreter','latex');
title('Interest Rate ($r$)', 'Fontsize', 15 ,'interpreter','latex');
%saveas(gcf,'InterestRate.png') % optionally: save plot 

% Suggest a case (1, 2 or 4)
fprintf('\nSuggested Case for the quarterly interest rate (r): Case 2\n')
fprintf(' -> We cannot detect a time trend (up- or downwards), so case 4 is eliminated.\n')
fprintf('    When deciding between case 1 and 2, case 2 is chosen because it seemes that the mean of the series seems to be unequal to zero.\n')


% plot quarterly CPI (second column of the data)
figure('Position', [0 0 650 450]) % fix size of the plot
plot(SWISS(:,2), "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('CPI', 'Fontsize', 14,'interpreter','latex');
title('Consumer Price Index ($p$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'CPI.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\nSuggested Case for the quarterly consumer price index (p): Case 4\n')
fprintf(' -> We see an upward-trending behaviour, so we decide on case 4.\n')


% plot quarterly GDP (third column)
figure('Position', [0 0 650 450]) % fix size of the plot
plot(SWISS(:,3), "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('Mio. CHF', 'Fontsize', 14, 'interpreter','latex');
title('Gross Domestic Product ($g$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'GDP.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\nSuggested Case for the quarterly GDP (g): Case 4\n')
fprintf(' -> We see an upward-trending behaviour, so we decide on case 4.\n')


% plot quarterly money stock (fourth column)
figure('Position', [0 0 650 450]) % fix size of the plot
plot(SWISS(:,4), "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('Mio. CHF', 'Fontsize', 14, 'interpreter','latex');
title('Money Stock M1 ($m$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'M1.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\nSuggested Case for the quarterly Money Stock (m): Case 4\n')
fprintf(' -> We see an upward-trending behaviour, so we decide on case 4.\n')

%% Exercise 1-3
fprintf('\n\n Exercise 1-3 (Section 3.2 in paper) ------------------------------------------------');

% -> the function is defined and explained in ADF.m (ADF test with k = 2)

% calling function for interest rate with case 2
[rho_hat_ir, se_rho_hat_ir, t_statistic_ir] = ADF(SWISS(:,1), "case2");

fprintf('\n\n------- Interest Rate -------')
fprintf('\nUsing case 2 (constant, no time trend in the estimated model) for')
fprintf('\nthe interest rate, we obtain the following estimates: \n')
fprintf('\nrho: %.4f\nStandard Error: %.4f\nt-statistic: %.4f\n',[rho_hat_ir, se_rho_hat_ir, t_statistic_ir]')



% calling function for CPI with case 4
[rho_hat_cpi, se_rho_hat_cpi, t_statistic_cpi] = ADF(SWISS(:,2), "case4");

fprintf('\n\n------- Consumer Price Index -------')
fprintf('\nUsing case 4 (constant, time trend in the estimated model) for')
fprintf('\nthe CPI, we obtain the following estimates: \n')
fprintf('\nrho: %.4f\nStandard Error: %.4f\nt-Statistic: %.4f\n',[rho_hat_cpi, se_rho_hat_cpi, t_statistic_cpi]')



% calling function for GDP with case 4
[rho_hat_gdp, se_rho_hat_gdp, t_statistic_gdp] = ADF(SWISS(:,3), "case4");

fprintf('\n\n------- Gross Domestic Product -------')
fprintf('\nUsing case 4 (constant, time trend in the estimated model) for')
fprintf('\nthe GDP, we obtain the following estimates: \n')
fprintf('\nrho: %.4f\nStandard Error: %.4f\nt-Statistic: %.4f\n',[rho_hat_gdp, se_rho_hat_gdp, t_statistic_gdp]')



% calling function for interest rate with case 2
[rho_hat_m1, se_rho_hat_m1, t_statistic_m1] = ADF(SWISS(:,4), "case4");

fprintf('\n\n------- Money Stock -------')
fprintf('\nUsing case 4 (constant, time trend in the estimated model) for')
fprintf('\nM1, we obtain the following estimates: \n')
fprintf('\nrho: %.4f\nStandard Error: %.4f\nt-Statistic: %.4f\n',[rho_hat_m1, se_rho_hat_m1, t_statistic_m1]')


%% Exercise 1-4

% -> the function is defined and explained in adv_ADF.m


%% Exercise 1-5
fprintf('\n\n Exercise 1-5 (Section 3.2 in paper) ------------------------------------------------');


%%%% Interest Rate with case 2 %%%%

% create empty matrix to store results
interest_rate_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(SWISS(:,1), "case2", i);
    
    % store the results
    interest_rate_results(i+1, 1) = k;
    interest_rate_results(i+1, 2) = rho;
    interest_rate_results(i+1, 3) = se_rho;
    interest_rate_results(i+1, 4) = t_statistic;
    interest_rate_results(i+1, 5) = SBC;
    interest_rate_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the interest rate:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(interest_rate_results);


%%%% CPI with case 4 %%%%

% create empty matrix to store results
CPI_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(SWISS(:,2), "case4", i);
    
    % store the results
    CPI_results(i+1, 1) = k;
    CPI_results(i+1, 2) = rho;
    CPI_results(i+1, 3) = se_rho;
    CPI_results(i+1, 4) = t_statistic;
    CPI_results(i+1, 5) = SBC;
    CPI_results(i+1, 6) = AIC;    
end

fprintf('\n\nResults depending on the lags (k) included for the CPI:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(CPI_results);



%%%% GDP with case 4 %%%%

% create empty matrix to store results
GDP_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(SWISS(:,3), "case4", i);
    
    % store the results
    GDP_results(i+1, 1) = k;
    GDP_results(i+1, 2) = rho;
    GDP_results(i+1, 3) = se_rho;
    GDP_results(i+1, 4) = t_statistic;
    GDP_results(i+1, 5) = SBC;
    GDP_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the GDP:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(GDP_results);




%%%% M1 with case 4 %%%%

% create empty matrix to store results
M1_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(SWISS(:,4), "case4", i);
    
    % store the results
    M1_results(i+1, 1) = k;
    M1_results(i+1, 2) = rho;
    M1_results(i+1, 3) = se_rho;
    M1_results(i+1, 4) = t_statistic;
    M1_results(i+1, 5) = SBC;
    M1_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the Money Stock:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(M1_results);


%% Exercise 1-6
fprintf('\n\n Exercise 1-6 (Section 3.2 in paper) ------------------------------------------------');


%%%% Interest Rate
% storing optimal results for interest rate based on SBC
opt_results_ir_SBC = interest_rate_results(interest_rate_results(:,5) == min(interest_rate_results(:,5)),:);

% storing optimal results for interest rate based on AIC
opt_results_ir_AIC = interest_rate_results(interest_rate_results(:,6) == min(interest_rate_results(:,6)),:);

fprintf('\n\n------- Interest Rate -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [opt_results_ir_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_ir_SBC(2), opt_results_ir_SBC(3), opt_results_ir_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [opt_results_ir_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_ir_AIC(2), opt_results_ir_AIC(3), opt_results_ir_AIC(4)]);



%%%% CPI
% storing optimal results for CPI based on SBC
opt_results_CPI_SBC = CPI_results(CPI_results(:,5) == min(CPI_results(:,5)),:);

% storing optimal results for interest rate based on AIC
opt_results_CPI_AIC = CPI_results(CPI_results(:,6) == min(CPI_results(:,6)),:);

fprintf('\n\n-------Consumer Price Index -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [opt_results_CPI_SBC(1)]); 
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_CPI_SBC(2), opt_results_CPI_SBC(3), opt_results_CPI_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [opt_results_CPI_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_CPI_AIC(2), opt_results_CPI_AIC(3), opt_results_CPI_AIC(4)]);



%%%% GDP
% storing optimal results for interest rate based on SBC
opt_results_GDP_SBC = GDP_results(GDP_results(:,5) == min(GDP_results(:,5)),:);

% storing optimal results for interest rate based on AIC
opt_results_GDP_AIC = GDP_results(GDP_results(:,6) == min(GDP_results(:,6)),:);

fprintf('\n\n------- Gross Domestic Product -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [opt_results_GDP_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_GDP_SBC(2), opt_results_GDP_SBC(3), opt_results_GDP_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [opt_results_GDP_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_GDP_AIC(2), opt_results_GDP_AIC(3), opt_results_GDP_AIC(4)]);



%%%% Money Stock
% storing optimal results for interest rate based on SBC
opt_results_M1_SBC = M1_results(M1_results(:,5) == min(M1_results(:,5)),:);

% storing optimal results for interest rate based on AIC
opt_results_M1_AIC = M1_results(M1_results(:,6) == min(M1_results(:,6)),:);

fprintf('\n\n------- Money Stock -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [opt_results_M1_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_M1_SBC(2), opt_results_M1_SBC(3), opt_results_M1_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [opt_results_M1_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [opt_results_M1_AIC(2), opt_results_M1_AIC(3), opt_results_M1_AIC(4)]);


fprintf('\n\nSBC and AIC do not always prefer the number of lags (k) because they penalize model \ncomplexity in a different way. ');
fprintf('In the present case with T = 156, SBC \npenalizes larger k stronger than AIC (for more details, see Section 3.2 in paper).');


%% Exercise 1-7
fprintf('\n\n Exercise 1.7 (Section 3.2 in paper) ------------------------------------------------');

fprintf('\n\nCritical Values of DF t-statistic at 5%% significance level: ');
fprintf('\nCase 1: -1.95');
fprintf('\nCase 2: -2.86');
fprintf('\nCase 4: -3.41');

fprintf('\n\n In the following, the SBC-preferred model specifications are considered: \n');

%%%% Interest Rate
fprintf('\n\n------- Interest Rate (Case 2) -------')
fprintf('\nt-Statistic: %.4f > -2.86', [opt_results_ir_SBC(4)]);
fprintf('\n -> We cannot reject the H0 of a unit root\n');

%%%% Consumer Prixe Index
fprintf('\n------- Consumer Prixe Index (Case 4) -------')
fprintf('\nt-Statistic: %.4f > -3.41', [opt_results_CPI_SBC(4)]);
fprintf('\n -> We cannot reject the H0 of a unit root\n');

%%%% Gross Domestic Product
fprintf('\n------- Gross Domestic Product (Case 4) -------')
fprintf('\nt-Statistic: %.4f > -3.41', [opt_results_GDP_SBC(4)]);
fprintf('\n -> We cannot reject the H0 of a unit root\n');

%%%% Money Stock
fprintf('\n------- Money Stock (Case 4) -------')
fprintf('\nt-Statistic: %.4f > -3.41', [opt_results_M1_SBC(4)]);
fprintf('\n -> We cannot reject the H0 of a unit root\n');

%% Exercise 1-8
fprintf('\n\n Exercise 1-8 (Section 3.3 in paper) ------------------------------------------------');

% transform the interest rate series using differences
delta_interest_rate = diff(SWISS(:,1));

% transform the CPI series using log differences
delta_CPI = diff(log(SWISS(:,2)));

% transform the GDP series using log differences
delta_GDP = diff(log(SWISS(:,3)));

% transform the Money M1 series using log differences
delta_M1 = diff(log(SWISS(:,4)));

fprintf('\n\nSeries successfully transformed using (log) differences. \n');

%% Exercise 1-9
fprintf('\n\n Exercise 1-9 (Section 3.3 in paper) ------------------------------------------------');

% plot delta quarterly interest rate
figure('Position', [0 0 650 450]) % fix size of the plot
plot(delta_interest_rate, "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('delta \%', 'Fontsize', 14,'interpreter','latex');
title('Tranformed Interest Rate ($\tilde{r}$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'DeltaInterestRate.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\n\nSuggested Case for the quarterly interest rate: Case 2\n')
fprintf(' -> We cannot detect a time trend (up- or downwards), so case 4 is eliminated.\n')
fprintf('    We see a reverting behaviour, so we decide for case 2.\n')


% plot delta quarterly CPI
figure('Position', [0 0 650 450]) % fix size of the plot
plot(delta_CPI, "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('delta ln CPI', 'Fontsize', 14,'interpreter','latex');
title('Transformed Consumer Price Index ($\tilde{p}$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'DeltaCPI.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\n\nSuggested Case for the quarterly Consumer Price Index: Case 2\n')
fprintf(' -> We cannot detect a time trend (up- or downwards), so case 4 is eliminated.\n')
fprintf('    We see a reverting behaviour, so we decide for case 2.\n')


% plot delta quarterly GDP
figure('Position', [0 0 650 450]) % fix size of the plot
plot(delta_GDP, "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('delta ln Mio. CHF', 'Fontsize', 14, 'interpreter','latex');
title('Transformed Gross Domestic Product ($\tilde{g}$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'DeltaGDP.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\n\nSuggested Case for the quarterly Gross Domestic Product: Case 2\n')
fprintf(' -> We cannot detect a time trend (up- or downwards), so case 4 is eliminated.\n')
fprintf('    We see a reverting behaviour, so we decide for case 2.\n')


% plot delta quarterly Money Stock
figure('Position', [0 0 650 450]) % fix size of the plot
plot(delta_M1, "b", 'linewidth', 1);
xlabel('t', 'Fontsize', 14,'interpreter','latex');
ylabel('delta ln Mio. CHF', 'Fontsize', 14, 'interpreter','latex');
title('Tranformed Money Stock M1 ($\tilde{m}$)', 'Fontsize', 15,'interpreter','latex');
%saveas(gcf,'DeltaM1.png') % optionally: save plot

% Suggest a case (1, 2 or 4)
fprintf('\n\nSuggested Case for the quarterly Money Stock: Case 2\n')
fprintf(' -> We cannot detect a time trend (up- or downwards), so case 4 is eliminated.\n')
fprintf('    We see a reverting behaviour, so we decide for case 2.\n')


%% Exercise 1-10 (Part 1-5)
fprintf('\n\n Exercise 1-10 (Part 1-5) (Section 3.3 in paper) ------------------------------------------------');


%%%% Transformed Interest Rate with case 2 %%%%

% create empty matrix to store results
delta_interest_rate_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(delta_interest_rate, "case2", i);
    
    % store the results
    delta_interest_rate_results(i+1, 1) = k;
    delta_interest_rate_results(i+1, 2) = rho;
    delta_interest_rate_results(i+1, 3) = se_rho;
    delta_interest_rate_results(i+1, 4) = t_statistic;
    delta_interest_rate_results(i+1, 5) = SBC;
    delta_interest_rate_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the transformed interest rate:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(delta_interest_rate_results);


%%%% Transformed CPI with case 2 %%%%

% create empty matrix to store results
delta_CPI_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(delta_CPI, "case2", i);
    
    % store the results
    delta_CPI_results(i+1, 1) = k;
    delta_CPI_results(i+1, 2) = rho;
    delta_CPI_results(i+1, 3) = se_rho;
    delta_CPI_results(i+1, 4) = t_statistic;
    delta_CPI_results(i+1, 5) = SBC;
    delta_CPI_results(i+1, 6) = AIC;    
end

fprintf('\n\nResults depending on the lags (k) included for the transformed CPI:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(delta_CPI_results);


%%%% Transformed GDP with case 2 %%%%

% create empty matrix to store results
delta_GDP_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(delta_GDP, "case2", i);
    
    % store the results
    delta_GDP_results(i+1, 1) = k;
    delta_GDP_results(i+1, 2) = rho;
    delta_GDP_results(i+1, 3) = se_rho;
    delta_GDP_results(i+1, 4) = t_statistic;
    delta_GDP_results(i+1, 5) = SBC;
    delta_GDP_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the transformed GDP:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(delta_GDP_results);


%%%% Transformed M1 with case 2 %%%%

% create empty matrix to store results
delta_M1_results = zeros(17, 6);

for i = 0:16 % iterate through k = 0, 1, ..., 16
    
    % calling the funtion depending on the lags
    [k, rho, se_rho, t_statistic, SBC, AIC] = adv_ADF(delta_M1, "case2", i);
    
    % store the results
    delta_M1_results(i+1, 1) = k;
    delta_M1_results(i+1, 2) = rho;
    delta_M1_results(i+1, 3) = se_rho;
    delta_M1_results(i+1, 4) = t_statistic;
    delta_M1_results(i+1, 5) = SBC;
    delta_M1_results(i+1, 6) = AIC;
end

fprintf('\n\nResults depending on the lags (k) included for the transformed M1:\n');
fprintf('\n       k       rho    Std. Error  t-stat.    SBC       AIC  \n');
disp(delta_M1_results);


%% Exercise 1-10 (Part 1-6)
fprintf('\n\n Exercise 1-10 (Part 1-6) (Section 3.3 in paper) ------------------------------------------------');


%%%% Interest Rate
% storing optimal results for delta interest rate based on SBC
delta_opt_results_ir_SBC = delta_interest_rate_results(delta_interest_rate_results(:,5) == min(delta_interest_rate_results(:,5)),:);

% storing optimal results for delta interest rate based on AIC
delta_opt_results_ir_AIC = delta_interest_rate_results(delta_interest_rate_results(:,6) == min(delta_interest_rate_results(:,6)),:);

fprintf('\n\n------- Interest Rate -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [delta_opt_results_ir_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_ir_SBC(2), delta_opt_results_ir_SBC(3), delta_opt_results_ir_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [delta_opt_results_ir_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_ir_AIC(2), delta_opt_results_ir_AIC(3), delta_opt_results_ir_AIC(4)]);



%%%% CPI
% storing optimal results for delta CPI based on SBC
delta_opt_results_CPI_SBC = delta_CPI_results(delta_CPI_results(:,5) == min(delta_CPI_results(:,5)),:);

% storing optimal results for delta CPI based on AIC
delta_opt_results_CPI_AIC = delta_CPI_results(delta_CPI_results(:,6) == min(delta_CPI_results(:,6)),:);

fprintf('\n\n-------Consumer Price Index -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [delta_opt_results_CPI_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_CPI_SBC(2), delta_opt_results_CPI_SBC(3), delta_opt_results_CPI_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [delta_opt_results_CPI_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_CPI_AIC(2), delta_opt_results_CPI_AIC(3), delta_opt_results_CPI_AIC(4)]);



%%%% GDP
% storing optimal results for delta GDP based on SBC
delta_opt_results_GDP_SBC = delta_GDP_results(delta_GDP_results(:,5) == min(delta_GDP_results(:,5)),:);

% storing optimal results for delta GDP based on AIC
delta_opt_results_GDP_AIC = delta_GDP_results(delta_GDP_results(:,6) == min(delta_GDP_results(:,6)),:);

fprintf('\n\n------- Gross Domestic Product -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [delta_opt_results_GDP_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_GDP_SBC(2), delta_opt_results_GDP_SBC(3), delta_opt_results_GDP_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [delta_opt_results_GDP_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_GDP_AIC(2), delta_opt_results_GDP_AIC(3), delta_opt_results_GDP_AIC(4)]);



%%%% Money Stock
% storing optimal results for delta M1 based on SBC
delta_opt_results_M1_SBC = delta_M1_results(delta_M1_results(:,5) == min(delta_M1_results(:,5)),:);

% storing optimal results for delta M1 based on AIC
delta_opt_results_M1_AIC = delta_M1_results(delta_M1_results(:,6) == min(delta_M1_results(:,6)),:);

fprintf('\n\n------- Money Stock -------')
fprintf('\nOptimal number of lags depending on SBC: %.0f', [delta_opt_results_M1_SBC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_M1_SBC(2), delta_opt_results_M1_SBC(3), delta_opt_results_M1_SBC(4)]);
fprintf('\nOptimal number of lags depending on AIC: %.0f', [delta_opt_results_M1_AIC(1)]);
fprintf('\nCorresponding Rho: %.4f, Standard Error %.4f, t-Statistic: %.4f\n', [delta_opt_results_M1_AIC(2), delta_opt_results_M1_AIC(3), delta_opt_results_M1_AIC(4)]);



%% Exercise 1-10 (Part 1-7)
fprintf('\n\n Exercise 1-10 (Part 1-7) (Section 3.3 in paper) ------------------------------------------------');

fprintf('\n\nCritical Values of DF t-statistic at 5%% significance level: ');
fprintf('\nCase 1: -1.95');
fprintf('\nCase 2: -2.86');
fprintf('\nCase 4: -3.41');


%%%% Interest Rate
fprintf('\n\n------- Interest Rate (Case 2) -------')
fprintf('\nt-Statistic: %.4f < -2.86', [delta_opt_results_ir_SBC(4)]);
fprintf('\n -> We can reject the H0 of a unit root\n');


%%%% Consumer Prixe Index
fprintf('\n------- Consumer Prixe Index (Case 2) -------')
fprintf('\nt-Statistic: %.4f > -3.41', [delta_opt_results_CPI_SBC(4)]);
fprintf('\n -> We cannot reject the H0 of a unit root\n');


%%%% Gross Domestic Product
fprintf('\n------- Gross Domestic Product (Case 2) -------')
fprintf('\nt-Statistic: %.4f < -3.41', [delta_opt_results_GDP_SBC(4)]);
fprintf('\n -> We can reject the H0 of a unit root\n');


%%%% Money Stock
fprintf('\n------- Money Stock (Case 2) -------')
fprintf('\nt-Statistic: %.4f < -3.41', [delta_opt_results_M1_SBC(4)]);
fprintf('\n -> We can reject the H0 of a unit root\n');



%% Exercise 1-11
fprintf('\n\n Exercise 1-11 (Section 3.3 in paper) ------------------------------------------------');

% Transformed GDP:
% obtain rho
rho_delta_GDP = delta_opt_results_GDP_SBC(2);

% obtain standard error of rho
se_rho_delta_GDP = delta_opt_results_GDP_SBC(3);

% compute upper and lower bound of CI with quantile of standard normal
% distribution (such that 95% CI)
CI_rho_GDP_lower = rho_delta_GDP - norminv(0.975)*se_rho_delta_GDP;
CI_rho_GDP_upper = rho_delta_GDP + norminv(0.975)*se_rho_delta_GDP;

fprintf('\n\n------- Confidence Intervall GDP -------')
fprintf('\n[%.4f, %.4f]\n', [CI_rho_GDP_lower, CI_rho_GDP_upper]);


% Transformed Money Stock M1:
% obtain rho
rho_delta_M1 = delta_opt_results_M1_SBC(2);

% obtain standard error of rho
se_rho_delta_M1 = delta_opt_results_M1_SBC(3);

% compute upper and lower bound of CI with quantile of standard normal
% distribution (such that 95% CI)
CI_rho_M1_lower = rho_delta_M1 - norminv(0.975)*se_rho_delta_M1;
CI_rho_M1_upper = rho_delta_M1 + norminv(0.975)*se_rho_delta_M1;

fprintf('\n\n------- Confidence Intervall Money Stock -------')
fprintf('\n[%.4f, %.4f]\n', [CI_rho_M1_lower, CI_rho_M1_upper]);


fprintf('\n The Confidence Intervalls displays all values for which we would not reject \n');
fprintf(' the null hypothesis at the given significane level (5%%).\n');


%% Exercise 2

% Exercises 2.1 - 2.3 are explained in the paper in section 4.1

%% Exercise 2-4

% include the VAR estimation toolbox
addpath("VAR2_toolbox");


%% Exercise 2-5

% set up the new dataset according to the ordering given by equation 3,
% i.e. Money Stock, Interest Rate, CPI, GDP
data = [delta_M1, delta_interest_rate, delta_CPI, delta_GDP];

% estimate the model including one lag
[constant, Phi_1, ~, Omega] = PHI(data,1);

%% Exercise 2-6
fprintf('\n\n Exercise 2-6 (Section 4.2 in paper) ------------------------------------------------\n');

% show the vector of estimates for the constant (c)
fprintf("\nestimate of c = \n")
disp(constant)

% show the estimator for the phi matrix (here phi_1 since we only inlcude
% one lag)
fprintf("\nestimate of phi_1 = \n")
disp(Phi_1)

%% Exercise 3

%% Exercise 3-1
fprintf('\n\n Exercise 3-1 (Section 4.2 in paper) ------------------------------------------------\n');

% show the Omega matrix (variance-covariance matrix of composite shocks)
fprintf("\nestimate of Omega = \n")
disp(Omega)

%% Exercise 3-2

% We are facing 16 unknown parameters but only 10 unique equations to solve
% the system -> underidentification of 6
% The Cholesky decomposition imposes 6 restrictions and hence, identifies
% the equation.
% For more detailed explanations, plese see the paper (section 4.2).


%% Exercise 3-3 
fprintf('\n\n Exercise 3-3 (Section 4.2 in paper) ------------------------------------------------\n');

% compute P from Omega using the Cholesky decomposition 
P = chol(Omega)';

% verify that the decomposition worked as it should
if sum(sum(P*P' - Omega)) - 1 + 1 == 0
    fprintf("\nDecomposition worked! \n")
else 
    fprintf("\nUps! Something didn't quite work out here... \n")
end

% show the P matrix
fprintf("\nP = \n")
disp(P)

%% Exercise 3-4 
fprintf('\n\n Exercise 3-4 (Section 4.2 in paper) ------------------------------------------------\n');

% obtain C as the diagonal matrix of the P matrix
C =  diag(diag(P));
fprintf('\nC =\n')
disp(C)

% obtain A using A = P*inv(sqrt(D))
A = P*inv(C);
fprintf('\nA = \n')
disp(A)

% obtain D using D = C*C'
D = C*C';
fprintf('\nestimate of D =\n')
disp(D)

% obtain B0 using using inv(B0) = A => B0 = inv(A)
B0 = inv(A);
fprintf('\nestimate of B_0 =\n')
disp(B0)
