clear
clc
close all

fprintf("Exercise 1 \n")
fprintf("------------------------------- \n\n")


%% Exercise 1-1
zeta1 = [0.2 0.1 0.1; 0.1 0.5 0; 0.1, -0.2, 0.2];
B = [0.2; 0; -0.8];
A = [1; -1; -1];
omega = [0.01 0 0.0001; 0 0.02 -0.0001; 0.0001, -0.0001, 0.002];

%% Exercise 1-2

phi1 = -B*A' + eye(3) + zeta1;
phi2 = -zeta1;

%% Exercise 1-3
phi1_pol = eye(3) - phi1 - phi2;

%% Exercise 1-4

fprintf('Rank of Phi(1): %d \n',rank(phi1_pol));

% hence, not full rank matrix (n = 3) -> h = 1

%% Exercise 1-5
% see VAR2.m

%% Exercise 1-6
rng(6);

y_t1 = VAR2(phi1, phi2, omega, 3000);
y_t2 = VAR2(phi1, phi2, omega, 3000);
y_t3 = VAR2(phi1, phi2, omega, 3000);
y_t4 = VAR2(phi1, phi2, omega, 3000);

T=1:3000;

figure

subplot(2,2,1)
plot(T,y_t1(:,1), T,y_t1(:,2), T,y_t1(:,3));
title('Iteration 1');
xlabel('t');
ylabel('y_t');
legend('home market', 'exchange rate', 'foreign market', 'location', 'best');

subplot(2,2,2)
plot(T,y_t2(:,1), T,y_t2(:,2), T,y_t2(:,3));
title('Iteration 2');
xlabel('t');
ylabel('y_t');

subplot(2,2,3)
plot(T,y_t3(:,1), T,y_t3(:,2), T,y_t3(:,3));
title('Iteration 3');
xlabel('t');
ylabel('y_t');

subplot(2,2,4)
plot(T,y_t4(:,1), T,y_t4(:,2), T,y_t4(:,3));
title('Iteration 4');
xlabel('t');
ylabel('y_t');

% How many stochastic trends exist?
% - exchange rate & foreign market are mirrored at home market
% - g = n - h = 2 stochastic trends exist 

%% Exercise 1-7
% see cv.m

%% Exercise 1-8
rng(6);

mat_cv = zeros(500,3);

for i = 1:500
    y_t = VAR2(phi1, phi2, omega, 3000); % data simulation
    mat_cv(i,:) = cv(y_t); % data estimation
end

% first column of mat_cv: constant

%% Exercise 1-9

[f2,xi2] = ksdensity(mat_cv(:,2)); 
[f3,xi3] = ksdensity(mat_cv(:,3)); 

figure
plot(xi2,f2, xi3,f3);
title('Visualization of the 2nd and 3rd element of the cv');
legend('Second', 'Third');

% Interpretation:
% - both densities very steep at 1 -> cv = [1 -1 -1]

%% Exercise 1-10
rng(6);
S = 2000;

ser1 = zeros(S,3);
eps = zeros(S,3);
eps(1,1) = 1;
for i = 1:S
    if i == 1
        ser1(i,:) = eps(i,:);
    elseif i == 2
        ser1(i,:) = phi1*ser1(i-1,:)' + eps(i,:)';
    else
        ser1(i,:) = phi1*ser1(i-1,:)' + phi2*ser1(i-2,:)' + eps(i,:)';
    end
end
col1 = ser1(S,:);

ser2 = zeros(S,3);
eps = zeros(S,3);
eps(1,2) = 1;
for i = 1:S
    if i == 1
        ser2(i,:) = eps(i,:);
    elseif i == 2
        ser2(i,:) = phi1*ser2(i-1,:)' + eps(i,:)';
    else
        ser2(i,:) = phi1*ser2(i-1,:)' + phi2*ser2(i-2,:)' + eps(i,:)';
    end
end
col2 = ser2(S,:);

ser3 = zeros(S,3);
eps = zeros(S,3);
eps(1,3) = 1;
for i = 1:S
    if i == 1
        ser3(i,:) = eps(i,:);
    elseif i == 2
        ser3(i,:) = phi1*ser3(i-1,:)' + eps(i,:)';
    else
        ser3(i,:) = phi1*ser3(i-1,:)' + phi2*ser3(i-2,:)' + eps(i,:)';
    end
end
col3 = ser3(S,:);

psi1_pol = [col1' col2' col3']

% Interpretation:
% - permanent impact of composite shocks on system variables
% - psi_21: permanent impact of u1 on y2
% - first: home-market price, second: exchange rate, third: foreign-market price

% - u2 has negative permanent impact on foreign market
% - u2 has strong positive permanent impact on exchange rate

%% Exercise 1-11

A_ort = null(A');
B_ort = null(B');

psi1_pol2 = A_ort * inv(B_ort'*(eye(3)-zeta1)*A_ort) * B_ort'
% same psi1 matrix as before with simulation

fprintf('Rank of Psi(1): %d \n\n\n', rank(psi1_pol));
% Intepretation:
% - reasonable as its rank is n-h = 3-1 = 2

%% Exercise 1-12

disp('Phi(1)*Psi(1) = ')
round(phi1_pol * psi1_pol, 10) % rounded to 10 digits

% Intepretation:
% - = 0 -> implication of cointegration on VAR(1)
% -> cointegration & computations work
% (in stationary setup: Phi(1)*Psi(1) = identity matrix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("\nExercise 2 \n")
fprintf("------------------------------- \n\n")

%% Exercise 2-1
% - need n*(n-1)/2 = 3 restrictions -> given by setting 3 elements of W = 0

% - contemporaneous effect of foreign market (3rd SV) on home market (1st SV) is ruled out
% - contemporaneous effect of home market (1rd SV) on exchange rate (2st SV) is ruled out
% - contemporaneous effect of foreign market (3rd SV) on exchange rate (12st SV) is ruled out

%% Exercise 2-2

% before: h, x, f -> now: x, h, f
omega
omega_r = zeros(3,3);
omega_r(1,1) = omega(2,2);
omega_r(1,2) = omega(1,2);
omega_r(1,3) = omega(2,3);
omega_r(2,1) = omega(2,1);
omega_r(2,2) = omega(1,1);
omega_r(2,3) = omega(1,3);
omega_r(3,1) = omega(3,2);
omega_r(3,2) = omega(3,1);
omega_r(3,3) = omega(3,3);
omega_r

%% Exercise 2-3
W_r = chol(omega_r, 'lower')
% or via: W_r = chol(omega_r)'

%% Exercise 2-4

W = zeros(3,3);
W(1,1) = W_r(2,2);
W(1,2) = W_r(2,1);
W(2,2) = W_r(1,1);
W(3,1) = W_r(3,2);
W(3,2) = W_r(3,1);
W(3,3) = W_r(3,3);

W

%% Exercise 2-5

rng(6);
S = 11;

ser1 = zeros(S,3);
u = zeros(S,3);
u(1,1) = 1;
for i = 1:S
    if i == 1
        ser1(i,:) = (W*u(i,:)')';
    elseif i == 2
        ser1(i,:) = phi1*ser1(i-1,:)' + W*u(i,:)';
    else
        ser1(i,:) = phi1*ser1(i-1,:)' + phi2*ser1(i-2,:)' + W*u(i,:)';
    end
end

ser2 = zeros(S,3);
u = zeros(S,3);
u(1,2) = 1;
for i = 1:S
    if i == 1
        ser2(i,:) = (W*u(i,:)')';
    elseif i == 2
        ser2(i,:) = phi1*ser2(i-1,:)' + W*u(i,:)';
    else
        ser2(i,:) = phi1*ser2(i-1,:)' + phi2*ser2(i-2,:)' + W*u(i,:)';
    end
end

ser3 = zeros(S,3);
u = zeros(S,3);
u(1,3) = 1;
for i = 1:S
    if i == 1
        ser3(i,:) = (W*u(i,:)')';
    elseif i == 2
        ser3(i,:) = phi1*ser3(i-1,:)' + W*u(i,:)';
    else
        ser3(i,:) = phi1*ser3(i-1,:)' + phi2*ser3(i-2,:)' + W*u(i,:)';
    end
end


x = 0:1:10;

figure

subplot(2,2,1)
plot(x,ser1(:,1),x,ser2(:,1),x,ser3(:,1))
title('Impact on home market (h)');
xlabel('s');
legend('Unit idiosyncratic shock in h (u1)','Unit idiosyncratic shock in m (u2)', 'Unit idiosyncratic shock in f (u3)');

subplot(2,2,2)
plot(x,ser1(:,2),x,ser2(:,2),x,ser3(:,2))
title('Impact on exchange rate (x)');
xlabel('s');

subplot(2,2,3)
plot(x,ser1(:,3),x,ser2(:,3),x,ser3(:,3))
title('Impact on foreign market (f)');
xlabel('s');

% Interpretation
% - permanent impacts of orthogonal shocks

