%{
matrix = [500 384 50 ;
          0   116 125 ;
          0   0   325];

close all
heat_map = heatmap(matrix);
heat_map.Title = 'how many times is a k present';
heat_map.XLabel = 'max K';
heat_map.YLabel = 'individual k';
%}
matrix = zeros(10);
column = 0;
for i = 1:10
    row = 0;
   column = column+1;
   if i == 1
       matrix(1,column) = 83300;
   else
       filename =  strcat('K_paths_values_', num2str(i), '_paths.mat');
       load(filename)
       for j = 1:10
            row = row+1;
            element = sum(sum(eval(strcat('all_routes_selected_', num2str(i))) == j));
            matrix(row,column) = element;
       end
   end
end
matrix = round(matrix./833,2);
x = [1 10];
y = [1 10];

close all

heat_map = heatmap(matrix);
heat_map.Title = 'Percentage of time the k-th path is optimal';
heat_map.XLabel = 'Maximum K';
heat_map.YLabel = 'Individual k';
heat_map.ColorScaling = 'log';
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.TickLabels = {'10','100'};

%{
close all
load('K_paths_values_10_paths.mat')
 
plot(1000:1000:140000,cost_reductions)
grid on

ax = ancestor(fig, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.0f');
ax.YAxis.Exponent = 0;
ytickformat('%.0f');
ax.ZAxis.Exponent = 0;
ztickformat('%.0f');
xticks(0:20000:140000)
title('Initial cost and optimized cost of upgrading for K = 10 and alpha = 0')
xlabel('Budget in millions of €')
ylabel('Reduction in millions of €')
legend('Cost reduction', 'location', 'northwest');
%}