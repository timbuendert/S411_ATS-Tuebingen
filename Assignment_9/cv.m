function [b_est] = cv(series)
T = size(series, 1);
y1 = series(:,1); % home market
y2 = series(:, 2:3); % exchange rate & foreign market

X = [ones(T,1), y2];
b_est = inv(X'*X)*(X'*y1);
end