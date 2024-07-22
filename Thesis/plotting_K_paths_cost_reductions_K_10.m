clear
close all

load('K_paths_values_10_paths.mat')
fig = plot(1000:1000:140000,cost_reductions);
ax = ancestor(fig, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.0f');
grid on
xlabel('Budget in millions of €')
ylabel('Cost reduction in millions of €')
title('Cost reduction for K = 10 and alpha = 0')
legend('Cost reduction')