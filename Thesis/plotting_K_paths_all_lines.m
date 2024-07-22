clear
budgets = 1000:1000:140000;

load('K_paths_values_1_path.mat')
load('K_paths_values_2_paths.mat')
load('K_paths_values_3_paths.mat')
load('K_paths_values_4_paths.mat')
load('K_paths_values_5_paths.mat')
load('K_paths_values_6_paths.mat')
load('K_paths_values_7_paths.mat')
load('K_paths_values_8_paths.mat')
load('K_paths_values_9_paths.mat')
load('K_paths_values_10_paths.mat')


close all
fig = plot(budgets,optimal_values);
hold on
plot(budgets,optimal_values_2)
hold on
plot(budgets,optimal_values_3)
hold on
plot(budgets,optimal_values_4)
hold on
plot(budgets,optimal_values_5)
hold on
plot(budgets,optimal_values_6)
hold on
plot(budgets,optimal_values_7)
hold on
plot(budgets,optimal_values_8)
hold on
plot(budgets,optimal_values_9)
hold on
plot(budgets,optimal_values_10)

ax = ancestor(fig, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.0f');
xticks(0:20000:140000)
title('Maximum shortest path per budget')
xlabel('Budget in millions of €')
ylabel('Maximum shortest path in hours')
grid on
 
legend('K = 1', 'K = 2', 'K = 3', 'K = 4', 'K = 5', 'K = 6', 'K = 7', 'K = 8', 'K = 9', 'K = 10')