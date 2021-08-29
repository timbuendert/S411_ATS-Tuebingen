function [series] = VAR2(phi1, phi2, omega, T)

series = zeros(T, 3);

eps = mvnrnd(zeros(3,1),omega,T); 
% or see hint in instructions:
% eps = chol(omega)'* normrnd(0,1,[T,3]); ??

for i = 1:T
    if i == 1
        series(i,:) = eps(i,:);
    elseif i == 2
        series(i,:) = phi1*series(i-1,:)' + eps(i,:)';
    else
        series(i,:) = phi1*series(i-1,:)' + phi2*series(i-2,:)' + eps(i,:)';
    end
end
