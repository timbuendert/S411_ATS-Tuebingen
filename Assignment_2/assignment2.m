clear
clc
close all

rng(200); % set seed to ensure replicability


%% Exercise 1-1
T = 100; % set number of observations, originally 100
AR1 = func_AR1(5,0.7,17,T); % call function based on parameters

figure;
plot(1:T, AR1, 'k'); % plot realization of process
title('Realization of AR(1)');
xlabel('i', 'Fontsize', 14);
ylabel('y', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 2-1
% -> see compute_CML_l_AR1()

%% Exercise 2-2

%AR1 = vertcat(17,AR1);  % add y0 to AR1 realization -> not required, L starts at t = 2

theta_a = [5, 0.7, 1]; % specify parameters for process a)
L_value_a = CLL_value(theta_a, AR1); % retrieve CLL value for a) given parameters and data

theta_b = [-5, 0.1, 0.25]; % specify parameters for process b)
L_value_b = CLL_value(theta_b, AR1); % retrieve CLL value for b) given parameters and data

fprintf('CLL Value of a): %f \nCLL Value of b): %f\n', L_value_a, L_value_b);

% Interpretation:
% - b) provides worse conditional log likelihood (CLL) value than a)
% -> reasonable given a) represents the "true" parameters (theta_0)
% -> decide for a) based on ML estimation

% CLL values always < 0: 
% -> ln(< 1) = negative
% (as 1/sqrt(2*pi*sigma_2) < 1 and exp(< 0) < 1 -> entire term < 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 3-1
% Outputs: 
% x: local minimum of the function -> parameter estimates
% fval: value of function at local minimum x
% exitflag: describes exit condition (1: successful convergence, 0: max. number 
% of iterations reached before solution, -1: algorithm terminated by output function)
% output: structure with information about optimization proces

func = @(theta)-1*CLL_value(theta, AR1); % create handle: set theta random and fix data, *-1 for minimization (or use CLL_value_min)
% Alternative: LogLike = @(x)-1*sum(log((1/(sqrt(2.*pi*x(3))))) - (((AR1(2:length(AR1)) - x(1) - (x(2).*AR1(1:length(AR1)-1))).^2)/(2.*x(3))));


% a) estimation with fminunc for varying start values
[estimates_1a, fval_1a] = fminunc(func, theta_a);
[estimates_1b, fval_1b] = fminunc(func, theta_b);

% b) estimation with fminsearch for varying start values
[estimates_2a, fval_2a, exitflag_2a, output_2a] = fminsearch(func, theta_a);
[estimates_2b, fval_2b, exitflag_2b, output_2b] = fminsearch(func, theta_b); 

% c) estimation with patternsearch for varying start values
options = optimoptions('patternsearch', 'MaxIterations', 500000); % increase number of iterations

[estimates_3a, fval_3a] = patternsearch(func, theta_a, [],[],[],[],[],[],[], options);
[estimates_3b, fval_3b] = patternsearch(func, theta_b, [],[],[],[],[],[],[], options);


% print results of minimization functions
estimates_a = [5, 0.7, 1, L_value_a; estimates_1a, -1*fval_1a; estimates_2a, -1*fval_2a; estimates_3a, -1*fval_3a];
estimates_b = [5, 0.7, 1, L_value_a; estimates_1b, -1*fval_1b; estimates_2b, -1*fval_2b; estimates_3b, -1*fval_3b];

function_names = ["theta_0      "; "fminunc      "; "fminsearch   "; "patternsearch"]';

fprintf('\n\n');
disp('With starting values from a):');
disp('                   phi          c       sigma_2  CLL Value');
fprintf('%s %9.3f  %10.3f %10.3f %10.3f \n',[function_names', estimates_a(:,1), estimates_a(:,2), estimates_a(:,3), estimates_a(:,4)]')

fprintf('\n\n');
disp('With starting values from b):');
disp('                   phi          c       sigma_2');
fprintf('%s %9.3f  %10.3f %10.3f %10.3f \n',[function_names', estimates_b(:,1), estimates_b(:,2), estimates_b(:,3), estimates_b(:,4)]')


% Interpretation
% - different optimization functions result in approximately same estimates with starting values a) and b) 
%   as should be since using same data (realization)
%   -> patternsearch stands out with slightly worse estimates with starting values b)
% - all result in slightly better CLL values than with true parameters in 2)


%% Exercise 3-2

phi_3a = estimates_1a(2)-0.25:0.0001:estimates_1a(2)+0.25; % set appropriate range for phi including estimate (0.62) and true parameter value (0.7)
phi_val_3a = zeros(size(phi_3a, 1),1); % initiate vector for CLL values
ind = 1; % initiate index
for i = phi_3a
    % use estimates_1a for c and sigma_2: from 3-1, fminunc, start values from a)
    L_phi = CLL_value([estimates_1a(1), i, estimates_1a(3)], AR1);
    phi_val_3a(ind, 1) = L_phi; % save CLL value in vector
    ind = ind + 1; % update index
end 

figure;
plot(phi_3a, phi_val_3a, 'k'); % plot realization of process
title('Conditional log-likelihood value with parameter estimates depending on phi for a)');
xlabel('phi', 'Fontsize', 14);
ylabel('CLL value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
xline(estimates_1a(2), "--r"); % add estimated phi as vertical line
xline(0.7, "--g"); % add true phi as vertical line

% Interpretation
% - using estimates for c and sigma_2: maximum CLL value is at phi = roughly 0.717
%   -> ML estimate of phi maixmizes L (very small deviation from phi_0 because not true c and sigma_2 and randomness)
% - peaked log likelihood function -> high estimation precision
