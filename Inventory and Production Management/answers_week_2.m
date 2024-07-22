%% Clear workspace and old figures
clear
close all
%% Load data
% Use csvread skipping the header row.  Alternatively, you can use the
% "Import data" wizard from the Home tab.
X = csvread('Lab1-2.csv', 1);
% Extract yearly demand, value, and current order quantities
D = X(:,1);
v = X(:,2);
Q = X(:,3);
% Number of items
n = length(D);
%% Current situation
% Use .* and ./ to get element-wise operations.
N_0 = sum(D ./ Q);                      % number of replenishments
TACS_0 = sum(0.5 * Q .* v);             % total average cycle stock
%% Derive the exchange curve
% TACS * N = c with c
c = 0.5 * sum(sqrt(D .* v))^2
% Reasonable lower and upper bounds for N
N_bounds = [0.3*N_0, 2*N_0];
% Make a vector of 100 different N values and compute the corresponding TACS values
N = linspace(N_bounds(1), N_bounds(2), 100);
TACS = c ./ N;
%% Make the figure
figure;
plot(N, TACS);
% Alternatively: fplot(@(n) c./n, N_bounds)
% Change some figure properties
hold on                                 % we are going to add more to this plot
grid
xlabel('Replenishments per year');
ylabel('Total average cycle stock');
title('Exchange curve');
xlim([0, max(xlim)]);                   % Set lower limit to 0, keep upper 
%unchanged
ylim([0, max(ylim)]);
% Add the current and new points
scatter(N_0, TACS_0, 'filled');
%% Compute points on the exchange curve
N_new(1) = N_0;
TACS_new(1) = c / N_0;
N_new(2) = c / TACS_0;
TACS_new(2) = TACS_0;
% Compute corresponding ratio A/r
Ar = TACS_new ./ N_new;
% Compute relative improvements
N_improve = (N_new - N_0) / N_0 * 100;
TACS_improve = (TACS_new - TACS_0) / TACS_0 * 100;
% Add points to the figure
scatter(N_new, TACS_new, 'filled');
legend('Exchange curve', 'Current position', 'Improved positions');
%% Print results
fprintf('Current position\n');
fprintf('N = %10.2f, TACS = %10.2f\n', N_0, TACS_0);
fprintf('\nNew positions\n');
fprintf('N = %10.2f (%6.2f %%), TACS = %10.2f (%6.2f %%), A/r = %6.2f\n', ...
        [N_new; N_improve; TACS_new; TACS_improve; Ar]);