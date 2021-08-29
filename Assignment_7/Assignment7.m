clear
clc
close all

fprintf("Exercise 1 \n")
fprintf("------------------------------- \n\n")


%% Exercise 1-1

data = xlsread('SWISS_1976_2014.xlsx');
% series 1: quarterly interest rate in percent (first column)
% series 2: quarterly Consumer Price Index (second column)
% series 3: quarterly Gross Domestic Product in Mio CHF (third column)
% series 4: quarterly Money Stock M1 in Mio CHF (fourth column)

%% Exercise 1-2
t = 1:1:size(data,1);
r = data(:,1);
p = data(:,2);
y = data(:,3); 
m = data(:,4);

figure

subplot(2,2,1)
plot(t,r)
title('Interest Rate ($r$)','interpreter','latex');
xlabel('t');
ylabel('%');

subplot(2,2,2)
plot(t,p)
title('Consumer Price Index ($p$)','interpreter','latex');
xlabel('t');
ylabel('Index value');

subplot(2,2,3)
plot(t,y)
title('Gross Domestic Product ($g$)','interpreter','latex');
xlabel('t');
ylabel('Mio. CHF');

subplot(2,2,4)
plot(t,m)
title('Money Stock M1 ($m$)','interpreter','latex');
xlabel('t');
ylabel('Mio. CHF');

% stylized facts: p, g, m exhibit trending structure
% all may contain unit root

%% Exercise 1-3

r_tilde = diff(data(:,1));
p_tilde = diff(log(data(:,2)));
y_tilde = diff(log(data(:,3)));
m_tilde = diff(log(data(:,4)));

t_tilde = 1:1:size(r_tilde,1);

figure

subplot(2,2,1)
plot(t_tilde,r_tilde)
title('Transformed Interest Rate ($\tilde{r}$)','interpreter','latex');
xlabel('t');
ylabel('%');

subplot(2,2,2)
plot(t_tilde,p_tilde)
title('Transformed Consumer Price Index ($\tilde{p}$)','interpreter','latex');
xlabel('t');
ylabel('Index value');

subplot(2,2,3)
plot(t_tilde,y_tilde)
title('Transformed Gross Domestic Product ($\tilde{g}$)','interpreter','latex');
xlabel('t');
ylabel('Mio. CHF');

subplot(2,2,4)
plot(t_tilde,m_tilde)
title('Transformed Money Stock M1 ($\tilde{m}$)','interpreter','latex');
xlabel('t');
ylabel('Mio. CHF');

% stylized facts: now all series exhibit strong fluctuations 
% -> seem more stationary

% advantage of using transformed data: rather stationary series (as also 
% unit root process is difference-stationary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 2-1
fprintf("Exercise 2 \n")
fprintf("------------------------------- \n\n")

addpath('VAR2_toolbox');

%% Exercise 2-2

x_t = [m_tilde, r_tilde, p_tilde, y_tilde]; % re-order according to VAR
[c,Phi,~,Omega]=PHI(x_t,1);

%% Exercise 2-3

[D,A,B_0,C] = Cholesky_decomposition(Omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 3-1
fprintf("Exercise 3 \n")
fprintf("------------------------------- \n\n")

[F,Psi] = VAR_2(Phi,0,10,1); % include only one lag
% first n=4 rows contain the matrix of parameter estimates for first lag 
% (forward iteration) etc.

%% Exercise 3-2

irf_o = orthogonalized_response(Psi,A,C,10,1);

% first column: effect of shock in first variable on first variable. 
% -> second column refers to a shock in the second variable on the first variable
% -> (n+1)th column refers to the effect of shock in first variable on second one

%% Exercise 3-3
x = 0:1:10;


figure

subplot(2,2,1)
plot(x,irf_o(:,1),x,irf_o(:,2),x,irf_o(:,3),x,irf_o(:,4))
title('Impact on m');
xlabel('s');
legend('One SD-shock in m','One SD-shock in r', 'One SD-shock in p', 'One SD-shock in y');

subplot(2,2,2)
plot(x,irf_o(:,5),x,irf_o(:,6),x,irf_o(:,7),x,irf_o(:,8))
title('Impact on r');
xlabel('s');

subplot(2,2,3)
plot(x,irf_o(:,9),x,irf_o(:,10),x,irf_o(:,11),x,irf_o(:,12))
title('Impact on p');
xlabel('s');

subplot(2,2,4)
plot(x,irf_o(:,13),x,irf_o(:,14),x,irf_o(:,15),x,irf_o(:,16))
title('Impact on y');
xlabel('s');


%% Exercise 3-4

% Describe and compare patterns of impulse response functions by answering the following questions: 
% How pronounced (big/small) are the responses to shocks, i.e., the response of one variable 
% to its own shocks and the shocks in the other variables? 
% How persistent are these shocks?

% - see effect of Cholesky ordering (allowing specific contemporaneous
% effects)
% - effect on own shock for all variables rather strong positive
% - all effect quickly die out -> indicate stationarity 
% - strong impacts for y and p, weak for r and in particular m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 4-1
fprintf("\nExercise 4\n")
fprintf("------------------------------- \n\n")

variancedec = variance_deco(Psi,A,D,10);

% first column: fraction of the variance of variable 1 that is due to 1
% second column: fraction of the variance of variable 2 that is due to 1
% (n+1)th column: fraction of the variance of variable 1 that is due to 2

%% Exercise 4-2

x = 1:1:10;


figure

subplot(2,2,1)
plot(x,variancedec(:,1),x,variancedec(:,5),x,variancedec(:,9),x,variancedec(:,13))
title('Variance of m');
xlabel('s');
legend('Variance share of m','Variance share of r', 'Variance share of p', 'Variance share of y', 'location', 'best');

subplot(2,2,2)
plot(x,variancedec(:,2),x,variancedec(:,6),x,variancedec(:,10),x,variancedec(:,14))
title('Variance of r');
xlabel('s');

subplot(2,2,3)
plot(x,variancedec(:,3),x,variancedec(:,7),x,variancedec(:,11),x,variancedec(:,15))
title('Variance of p');
xlabel('s');

subplot(2,2,4)
plot(x,variancedec(:,4),x,variancedec(:,8),x,variancedec(:,12),x,variancedec(:,16))
title('Variance of y');
xlabel('s');


%% Exercise 4-3

% Discuss the results and draw a conclusion from the plots of the variance decompositions. 
% How do the variances of the four series decompose?

% for m: basicially all variance attributed to m
% for r: significantly most attributed to r, some to m and nothing to others
% for p: basicially all variance attributed to p, some to y
% for y: almost half attributed to y and p, only some to r and m
% -> m, r, and p mostly attribute variance to themselves
% -> variance of y decomposed in mainly 2 variables

%% Exercise 4-4

x_t_new = [p_tilde, y_tilde, m_tilde, r_tilde]; % re-order 

[c_new,Phi_new,~,Omega_new]=PHI(x_t_new,1);
[D_new,A_new,B_0_new,C_new] = Cholesky_decomposition(Omega_new);
[F_new,Psi_new] = VAR_2(Phi_new,0,10,1);

irf_o = orthogonalized_response(Psi_new,A_new,C_new,10,1);
x = 0:1:10;
figure

subplot(2,2,1)
% legend('One SD-shock in p','One SD-shock in y', 'One SD-shock in m', 'One SD-shock in r');
% plot(x,irf_o(:,9),x,irf_o(:,10),x,irf_o(:,11),x,irf_o(:,12))
plot(x,irf_o(:,11),x,irf_o(:,12),x,irf_o(:,9),x,irf_o(:,10))
legend('One SD-shock in m','One SD-shock in r', 'One SD-shock in p', 'One SD-shock in y'); % to achieve same ordering as above
title('Impact on m');
xlabel('s');

subplot(2,2,2)
% plot(x,irf_o(:,13),x,irf_o(:,14),x,irf_o(:,15),x,irf_o(:,16))
plot(x,irf_o(:,15),x,irf_o(:,16),x,irf_o(:,13),x,irf_o(:,14))
title('Impact on r');
xlabel('s');

subplot(2,2,3)
%plot(x,irf_o(:,1),x,irf_o(:,2),x,irf_o(:,3),x,irf_o(:,4))
plot(x,irf_o(:,3),x,irf_o(:,4),x,irf_o(:,1),x,irf_o(:,2))
title('Impact on p');
xlabel('s');

subplot(2,2,4)
% plot(x,irf_o(:,5),x,irf_o(:,6),x,irf_o(:,7),x,irf_o(:,8))
plot(x,irf_o(:,7),x,irf_o(:,8),x,irf_o(:,5),x,irf_o(:,6))
title('Impact on y');
xlabel('s');

variancedec = variance_deco(Psi,A,D,10);
x = 1:1:10;
figure

subplot(2,2,1)
plot(x,variancedec(:,11),x,variancedec(:,15),x,variancedec(:,3),x,variancedec(:,7))
title('Variance of m');
xlabel('s');
% legend('Variance share of p','Variance share of y', 'Variance share of m', 'Variance share of r', 'location', 'best');
legend('Variance share of m','Variance share of r', 'Variance share of p', 'Variance share of y', 'location', 'best'); % to achieve same ordering as above

subplot(2,2,2)
plot(x,variancedec(:,12),x,variancedec(:,16),x,variancedec(:,4),x,variancedec(:,8))
title('Variance of r');
xlabel('s');

subplot(2,2,3)
plot(x,variancedec(:,9),x,variancedec(:,13),x,variancedec(:,1),x,variancedec(:,5))
title('Variance of p');
xlabel('s');

subplot(2,2,4)
plot(x,variancedec(:,10),x,variancedec(:,14),x,variancedec(:,2),x,variancedec(:,6))
title('Variance of y');
xlabel('s');

% How does this change affect your impulse response function and the 
% variance decomposition? Briefly explain the reason for this result.

% IRF
% - different ordering has effect on contemporaenous effect (per construction)
% - apart from that: same pattern of responses, also die out

% Variance Decomposition
% - only slight changes, similar pattern -> more changes for y and r

% -> reason: different ordering influences contemporaneous effects, but
% apart from that not much -> ordering does not matter a lot in this
% application
% -> economic interpretation / reason: (??) due to autoregressive dynamics, ordering
% not important

%% Exercise 4-5

% try out all different permutations: if that leads to very different 
% results -> ordering matters 
% best out-of-sample of different orderings (??)

% Ordering does not matter -> contemporaneous effects 0 
% -> not the case comparing B_0 and B_0_new

