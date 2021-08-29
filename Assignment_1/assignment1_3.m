clear
clc
close all

rng(103) % set seed to ensure replicability


%% Exercise 1
% -> see func_AR2Fmat()
[F1, Eig1] = func_AR2_Fmat(0.85,0.3,-0.2); % call function based on parameters in a)
[F2, Eig2] = func_AR2_Fmat(1,-0.2,0.2); % call function based on parameters in b)
[F3, Eig3] = func_AR2_Fmat(0.4,-0.5,0.5); % call function based on parameters in c)


%% Exercise 2
x = 1:100; % create sequences of js

% a)
F1_11 = zeros(100,1); % initialize matrix for the Fj_11 values
for i = x % iterate through all j
    F1_j = F1^i; % take F to the jth power
    F1_11(i,1) = F1_j(1,1); % read out its first column, first row element
end

figure;
plot(x, F1_11); % plot j and Fj_11 sequences
title('F_j11 for AR(3) in part a)');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - MA coefficients converge to 0 as t increases


% b) -> see comments of 2.a)
F2_11 = zeros(100,1);
for i = x
    F2_j = F2^i;
    F2_11(i,1) = F2_j(1,1);
end

figure;
plot(x, F2_11);
title('F_j11 for AR(3) in part b)');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - MA coefficients converge to 0.833 as t increases -> impact of past never fades


% c) -> see comments of 2.a)
F3_11 = zeros(100,1);
for i = x
    F3_j = F3^i;
    F3_11(i,1) = F3_j(1,1);
end

figure;
plot(x, F3_11);
title('F_j11 for AR(3) in part c)');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);
% Interpretation
% - MA coefficients converge to 0 after fluctuating into positive and negative dimensions as t increases


%% Exercise 3
% a)
F1_11_2 = zeros(100,1); % initialize second matrix for the Fj_11 values
[F1_eigVec, F1_eigenVal] = eig(F1); % retrieve eigenvectors and eigenvalues of F
for i = x % % iterate through all j
    F1_j_2 = F1_eigVec * F1_eigenVal^i * (eye(3)/F1_eigVec); % compute Fj using eigenvectors and - values
    F1_11_2(i,1) = F1_j_2(1,1); % read out its first column, first row element
end

figure;
plot(x, F1_11, x, F1_11_2, 'r'); % plot the two differenct computations for Fj_11 for comparison
title('Comparison of F_j11 for AR(3) in part a)');
legend('F Mult.', 'Eig. Decom.', 'Location', 'Northwest');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% b) -> see comments of 3.a)
F2_11_2 = zeros(100,1);
[F2_eigVec, F2_eigenVal] = eig(F2); 
for i = x
    F2_j_2 = F2_eigVec * F2_eigenVal^i * (eye(3)/F2_eigVec);
    F2_11_2(i,1) = F2_j_2(1,1);
end
figure;
plot(x, F2_11, x, F2_11_2, 'r');
title('Comparison of F_j11 for AR(3) in part b)');
legend('F Mult.', 'Eig. Decom.', 'Location', 'Northwest');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% c) -> see comments of 3.a)
F3_11_2 = zeros(100,1);
[F3_eigVec, F3_eigenVal] = eig(F3); 
for i = x
    F3_j_2 = F3_eigVec * F3_eigenVal^i * (eye(3)/F3_eigVec);
    F3_11_2(i,1) = F3_j_2(1,1);
end
figure;
plot(x, F3_11, x, F3_11_2, 'r');
title('Comparison of F_j11 for AR(3) in part c)');
legend('F Mult.', 'Eig. Decom.', 'Location', 'Northeast');
xlabel('j', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


% Interpretation
% - exactly same results for Fj comparing F Mult. and Eig. Decom.


%% Exercise 4
% -> see func_lagPol3()


%% Exercise 5
% a)
Pol1_range = (-3:0.1:3); % set plausible range for the polynomial
Pol1_eval = func_lagPol3(0.85, 0.3, -0.2, Pol1_range); % call function for parameters of a)
figure;
plot(Pol1_range, Pol1_eval, Pol1_range, zeros(size(Pol1_range,2),1), 'k'); % plot the evaluations of the polynomial and a 0-line
title('Polynomial for AR(3) in part a)');
xlabel('z', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% b) -> see comments of 5.a)
Pol2_range = (-2:0.1:2);
Pol2_eval = func_lagPol3(1, -0.2, 0.2, Pol2_range);
figure;
plot(Pol2_range, Pol2_eval, Pol2_range, zeros(size(Pol2_range,2),1), 'k');
title('Polynomial for AR(3) in part b)');
xlabel('z', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);

% c) -> see comments of 5.a)
Pol3_range = (-2:0.1:2);
Pol3_eval = func_lagPol3(0.4, -0.5, 0.5, Pol3_range);
figure;
plot(Pol3_range, Pol3_eval, Pol3_range, zeros(size(Pol3_range,2),1), 'k');
title('Polynomial for AR(3) in part c)');
xlabel('z', 'Fontsize', 14);
ylabel('value', 'Fontsize', 14);
set(gca, 'Fontsize', 12);


%% Exercise 6

% eigenvalues of F = inverses of roots of lag polynomial

% a)
% all eigenvalues real and < 1 in absolute value
% MA coefficients drop to 0 as j -> infinity
% 3 real roots of lag polynomial > 1 in absolute value
% -> stationary

% b)
% 1 real eigenvalue with value 1 (unit-root) / 2 complex eigenvalues with R < 1
% MA coefficients take on constant value as j -> infinity
% 1 real root of lag polynomial = 1
% -> not stationary

% c)
% 1 real eigenvalue < 1 in absolute value / 2 complex eigenvalues with R < 1
% MA coefficients drop to 0 as j -> infinity
% 1 real root of lag polynomial > 1 in absolute value
% -> stationary