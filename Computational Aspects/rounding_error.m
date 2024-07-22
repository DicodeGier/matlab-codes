clear
close all

% Step sizes
n = 13;
h = 10.^-(0:n-1)';

% Calculate finite difference approximation and errors
dfx0 = cos(0.5);
d = (sin(0.5+h) - sin(0.5-h)) ./ (2*h);
e = abs(dfx0 - d);
[~, i] = min(e);

% Plot errors
figure;
loglog(h, e, 'b');
hold on
scatter(h(i), e(i), 50, 'b', 'fill');
xlabel('step size');
ylabel('error');
title('Plot of Error vs. Step Size');
grid

% Now do it for some other points as well
dfx0 = cos(1);
d(:,2) = (sin(1+h) - sin(1-h)) ./ (2*h);
e(:,2) = abs(dfx0 - d(:,2));
[~, i] = min(e(:,2));
loglog(h, e(:,2), 'g');
scatter(h(i), e(i,2), 50, 'g', 'fill');

dfx0 = cos(2);
d(:,3) = (sin(2+h) - sin(2-h)) ./ (2*h);
e(:,3) = abs(dfx0 - d(:,3));
[~, i] = min(e(:,3));
loglog(h, e(:,3), 'r');
scatter(h(i), e(i,3), 50, 'r', 'fill');

dfx0 = cos(-0.2);
d(:,4) = (sin(-0.2+h) - sin(-0.2-h)) ./ (2*h);
e(:,4) = abs(dfx0 - d(:,4));
[~, i] = min(e(:,4));
loglog(h, e(:,4), 'm');
scatter(h(i), e(i,4), 50, 'm', 'fill');

