clear
clc

%%%parameters that can be changed%%%
cost_km_shsr = 16.5;
cost_km_hsr = 25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%load data and construct parameters
data = importfile('data network complete.xlsx', 'tijden (uren)');
hours = data.data;
clear data

data = importfile('data network complete.xlsx', 'afstanden (kilometers)');
distances = data.data;
clear data

%check if data is complete
for i = 1:size(hours,1)
    for j = i+1:size(hours,1)
        if hours(i,j) > 0
            if distances(i,j) > 0
                continue
            else
                fprintf('error at %g, %g\n',i,j)
            end
        end
    end
end

indices = triu(ones(size(hours,1)),1);
indices = indices > 0;

data = importfile('data network complete.xlsx', 'u_shsr');
u_shsr = data.data;
u_shsr(isnan(u_shsr)) = 0;
clear data
u_shsr_matrix = u_shsr;
u_shsr(~indices) = Inf;
u_shsr = u_shsr';
u_shsr = u_shsr(:);
u_shsr(u_shsr == Inf) = [];
u_shsr = u_shsr';

data = importfile('data network complete.xlsx', 'u_hsr');
u_hsr = data.data;
u_hsr(isnan(u_hsr)) = 0;
clear data
u_hsr_matrix = u_hsr;
u_hsr(~indices) = Inf;
u_hsr = u_hsr';
u_hsr = u_hsr(:);
u_hsr(u_hsr == Inf) = [];
u_hsr = u_hsr';

distances_matrix = distances;
distances(~indices) = Inf;
distances = distances';
distances = distances(:);
distances(distances == Inf) = [];
distances = distances';

t = hours;
t(~indices) = Inf;
t = t';
t = t(:);
t(t == Inf) = [];
t = t';

shsr = u_shsr .* distances./250;
hsr = u_hsr .* distances./300;

cost_shsr = cost_km_shsr.*distances.*u_shsr;
cost_hsr = cost_km_hsr.*distances.*u_hsr;

u_shsr_0 = u_shsr == 0;
u_shsr_1 = u_shsr == 1;
u_shsr_copy = u_shsr;
u_shsr_copy(u_shsr_0) = 1;
u_shsr_copy(u_shsr_1) = 0;
shsr = shsr + u_shsr_copy .* t;

u_hsr_0 = u_hsr == 0;
u_hsr_1 = u_hsr == 1;
u_hsr_copy = u_hsr;
u_hsr_copy(u_hsr_0) = 1;
u_hsr_copy(u_hsr_1) = 0;
hsr = hsr + u_hsr_copy .* t;

shsr_matrix = zeros(size(hours,1));
shsr_copy = shsr;
for i = 1:size(hours,1)
   for j = 1:size(hours,1)
       if j <= i
           shsr_matrix(i,j) = 0;
       else
           shsr_matrix(i,j) = shsr_copy(1);
           shsr_copy(1) = [];
       end
   end
end

hsr_matrix = zeros(size(hours,1));
hsr_copy = hsr;
for i = 1:size(hours,1)
   for j = 1:size(hours,1)
       if j <= i
           hsr_matrix(i,j) = 0;
       else
           hsr_matrix(i,j) = hsr_copy(1);
           hsr_copy(1) = [];
       end
   end
end

red_shsr = t - shsr;
red_hsr = t - hsr;

red_shsr_matrix = zeros(size(hours,1));
red_shsr_copy = red_shsr;
for i = 1:size(hours,1)
   for j = 1:size(hours,1)
       if j <= i
           red_shsr_matrix(i,j) = 0;
       else
           red_shsr_matrix(i,j) = red_shsr_copy(1);
           red_shsr_copy(1) = [];
       end
   end
end
   
red_hsr_matrix = zeros(size(hours,1));
red_hsr_copy = red_hsr;
for i = 1:size(hours,1)
   for j = 1:size(hours,1)
       if j <= i
           red_hsr_matrix(i,j) = 0;
       else
           red_hsr_matrix(i,j) = red_hsr_copy(1);
           red_hsr_copy(1) = [];
       end
   end
end


u_hsr = u_hsr';
u_shsr = u_shsr';

%%%construct the shortest paths
nodes = 49;
real_nodes = 35; %nodes between which the shortest path is calculated (i.e. not the artificial nodes)
possibilities = real_nodes*(real_nodes - 1)/2; %maximum number of shortest paths
result = cell(possibilities,6);

counter = 1;
for i = 1:nodes
    %exclude artificial nodes
    if i == 5 || i == 6 || i == 7 || i == 9 || i == 16 || i == 17 || i == 18 || i == 23 || i == 24 || i == 29 || i == 30 || i == 31 || i == 32 || i == 44
        continue
    end
    for j = i+1:nodes
        if j == 5 || j == 6 || j == 7 || j == 9 || j == 16 || j == 17 || j == 18 || j == 23 || j == 24 || j == 29 || j == 30 || j == 31 || j == 32 || j == 44
            continue
        end
        %calculate shortest path
        [cost, route] = dijkstra(hours,i,j);
        
        %store required information in a cell
        route = flip(route,2);
        result{counter,1} = i; %starting node
        result{counter,2} = j; %ending node
        result{counter,3} = route; %traversed route
        result{counter,4} = cost; %cost of route
        matrix = zeros(nodes);
        for k = 1:length(route)-1
           element_1 = route(k);
           element_2 = route(k+1);
           row_node = min(element_1,element_2);
           column_node = max(element_1,element_2);
           matrix(row_node,column_node) = 1;
        end
        result{counter,5} = matrix; %matrix in which edges are in the upper triangular part
        matrix(~indices) = Inf;
        matrix = matrix';
        matrix = matrix(:);
        matrix(matrix==Inf) = [];
        matrix = matrix';
        result{counter,6} = matrix; %vector of edges stretched along the second dimension, i.e. x = [x11 x12 x13 ...] so "the rows first"
        counter = counter + 1;
        
    end
end

%hist([result{:,4}])
max_length = max([result{:,4}]);

%%%GEHARDCODED: DIT NOG DYNAMISCHER MAKEN!!!!%%%
%exclude routes that are less than the mean travel time
number_routes = size(result,1);
number_edges = length(result{1,6});
x_st = zeros(number_routes,number_edges);
for k = 1:number_routes
    x_st(k,:) = result{k,6};
end

%objective function
z_obj = 1;
y_1_obj = zeros(1,number_edges);
y_2_obj = zeros(1,number_edges);

obj = [z_obj y_1_obj y_2_obj];

%first constraint: shortest path longer than any path after upgrading
y_1_first_constraint = red_shsr.*x_st;
y_2_first_constraint = red_hsr.*x_st;
z_first_constraint = ones(number_routes,1);
first_constraint_matrix = [z_first_constraint y_1_first_constraint y_2_first_constraint];
rhs_first_constraint = t*transpose(x_st);

%second constraint: edge can only be upgraded once
y_1_second_constraint = eye(number_edges);
y_2_second_constraint = eye(number_edges);
z_second_constraint = zeros(number_edges,1);
second_constraint_matrix = [z_second_constraint y_1_second_constraint y_2_second_constraint];
rhs_second_constraint = ones(number_edges,1);

%third constraint: only upgradeable edges are upgraded to shsr
z_third_constraint = zeros(number_edges,1);
y_1_third_constraint = eye(number_edges);
y_2_third_constraint = zeros(number_edges);
third_constraint_matrix = [z_third_constraint y_1_third_constraint y_2_third_constraint];
rhs_third_constraint = u_shsr;

%fourth constraint: only upgradeable edges are upgraded to hsr
z_fourth_constraint = zeros(number_edges,1);
y_1_fourth_constraint = zeros(number_edges);
y_2_fourth_constraint = eye(number_edges);
fourth_constraint_matrix = [z_fourth_constraint y_1_fourth_constraint y_2_fourth_constraint];
rhs_fourth_constraint = u_hsr;


%fifth constraint: cost within budget
z_fifth_constraint = 0;
y_1_fifth_constraint = cost_shsr;
y_2_fifth_constraint = cost_hsr;
fifth_constraint_matrix = [z_fifth_constraint y_1_fifth_constraint y_2_fifth_constraint];

constraint_matrix = [first_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; fifth_constraint_matrix];
sparse_constraint_matrix = sparse(constraint_matrix);
rhs_till_now = [transpose(rhs_first_constraint); rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint];

alphas = 0:0.002:0.05;
budgets = 1000:1000:140000;
z = zeros(length(alphas),length(budgets));

total_cost_upgrades_to_shsr = zeros(nodes);
total_cost_upgrades_to_hsr = zeros(nodes);

column = 0;
for budget = 1000:1000:140000 
    column = column+1;
    row = 0;
    for alpha = 0:0.002:0.05
        row = row+1;
        [optimal_value,used_budget,new_budget,~,~,cost_upgrade_matrix_shsr,cost_upgrade_matrix_hsr] = one_shortest_path_shortened(budget,alpha,nodes,rhs_till_now,obj,sparse_constraint_matrix,number_routes,number_edges,cost_shsr,cost_hsr);
        z(row,column) = used_budget - new_budget;
        total_cost_upgrades_to_shsr = total_cost_upgrades_to_shsr + cost_upgrade_matrix_shsr;
        total_cost_upgrades_to_hsr = total_cost_upgrades_to_hsr + cost_upgrade_matrix_hsr;
    end
end
%}
%load('one_shortest_path_cost_reductions.mat')
z(1,78) = 0;
z = round(z,2);
close all

fig = surf(budgets,alphas,z);
ax = ancestor(fig, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.0f');
ax.ZAxis.Exponent = 0;
ztickformat('%.0f');
xticks(0:20000:140000)
yticks(0:0.01:0.05)
title('Reduction of costs per budget and alpha')
xlabel('Budget in millions of €')
ylabel('alpha')
zlabel('reduction in costs in millions of €')
hColorbar = colorbar;
set(hColorbar, 'YTickLabel', cellstr(num2str(reshape(get(hColorbar, 'YTick'),[],1),'%0.0f')) )
grid on


cities = cell(size(hours,1),1);
cities{1} = 'Amsterdam';
cities{2} = 'Berlijn';
cities{3} = 'Hamburg';
cities{4} = 'München';
cities{5} = 'Frankfurt';
cities{6} = 'Hannover';
cities{7} = 'Keulen';
cities{8} = 'Brussel';
cities{9} = 'Lille';
cities{10} = 'Kopenhagen';
cities{11} = 'Oslo';
cities{12} = 'Stockholm';
cities{13} = 'Parijs';
cities{14} = 'Lyon';
cities{15} = 'Marseille';
cities{16} = 'Straatsburg';
cities{17} = 'Bordeaux';
cities{18} = 'Perpignan';
cities{19} = 'Bern';
cities{20} = 'Barcelona';
cities{21} = 'Madrid';
cities{22} = 'Sevilla';
cities{23} = 'Burgos';
cities{24} = 'Cordoba';
cities{25} = 'Praag';
cities{26} = 'Milaan';
cities{27} = 'Rome';  
cities{28} = 'Napels';
cities{29} =    'Venetië';
cities{30} =    'Bologna';
cities{31} =    'Florence';
cities{32} =    'Linz';
cities{33} =    'Wenen';
cities{34} =    'Ljubljana';
cities{35} =    'Zagreb'; 
cities{36} =    'Boedapest'; 
cities{37} =    'Warschau'; 
cities{38} =    'Katowice';
cities{39} =    'Belgrado';
cities{40} =    'Podgorica';
cities{41} =    'Boekarest';
cities{42} =    'Sofia';
cities{43} =    'Skopje';
cities{44} =    'Thessaloniki';
 cities{45} =   'Athene';
 cities{46} =   'Istanbul';
 cities{47} =   'Londen';
 cities{48} =   'Birmingham';
 cities{49} =   'Glasgow';
 
 total_edges_shsr = sum(sum(total_cost_upgrades_to_shsr > 0));
 cell_edges_shsr = cell(total_edges_shsr,3);
 counter = 1;
 for i = 1:49
     for j = 1:49
     if total_cost_upgrades_to_shsr(i,j) > 0
        cell_edges_shsr{counter,1} = cities{i};
        cell_edges_shsr{counter,2} = cities{j};
        cell_edges_shsr{counter,3} = total_cost_upgrades_to_shsr(i,j);
        counter = counter + 1;
     end
     end
 end
 
 sorted_cell_edges_shsr = sortrows(cell_edges_shsr,3);
 
 total_edges_hsr = sum(sum(total_cost_upgrades_to_hsr > 0));
 cell_edges_hsr = cell(total_edges_hsr,3);
 counter = 1;
 for i = 1:49
     for j = 1:49
     if total_cost_upgrades_to_hsr(i,j) > 0
        cell_edges_hsr{counter,1} = cities{i};
        cell_edges_hsr{counter,2} = cities{j};
        cell_edges_hsr{counter,3} = total_cost_upgrades_to_hsr(i,j);
        counter = counter + 1;
     end
     end
 end
 
 sorted_cell_edges_hsr = sortrows(cell_edges_hsr,3);