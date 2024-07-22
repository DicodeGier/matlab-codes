clear
clc
addpath("C:\gurobi1001\win64\matlab")
 
%%%parameters that can be changed%%%
cost_km_shsr = 16.5;
cost_km_hsr = 25;
budget = 26000;
alpha = 0.001;
K = 30;
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

t_zero = hours;
t_zero(~indices) = Inf;
t_zero = t_zero';
t_zero = t_zero(:);
t_zero(t_zero == Inf) = [];
t_zero = t_zero';

hours(hours == 0) = Inf;
t = hours;
t(~indices) = 0;
t = t';
t = t(:);
t(t == 0) = [];
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
shsr = shsr + u_shsr_copy .* t_zero;

u_hsr_0 = u_hsr == 0;
u_hsr_1 = u_hsr == 1;
u_hsr_copy = u_hsr;
u_hsr_copy(u_hsr_0) = 1;
u_hsr_copy(u_hsr_1) = 0;
hsr = hsr + u_hsr_copy .* t_zero;

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

red_shsr = t_zero - shsr;
red_hsr = t_zero - hsr;

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

results = {};
all_K = [];
counter = 1 ;   
for i = 1:size(hours,2)
    disp(i)
    for j = i+1:size(hours,2)
        [paths, lengths] = kShortestPath(hours, i, j, K);
        k_i_j = length(lengths);
        all_K(1,end+1) = k_i_j;
        for k = 1:k_i_j
            results{counter,1} = i;
            results{counter,2} = j;
            path = cell2mat(paths(1,k));
            results{counter,3} = path;
            results{counter,4} = lengths(1,k);
            matrix = zeros(size(hours,1));
            for l = 1:length(path) - 1
                element_1 = path(l);
                element_2 = path(l+1);
                row_element = min([element_1, element_2]);
                column_element = max([element_1, element_2]);
                matrix(row_element, column_element) = 1;
            end
            results{counter,5} = matrix;
            matrix(~indices) = Inf;
            matrix = matrix';
            matrix = matrix(:);
            matrix(matrix==Inf) = [];
            matrix = matrix';
            results{counter,6} = matrix;
            counter = counter + 1;
        end
    end
end

number_routes = size(results,1);
number_edges = length(results{1,6});
x_st = zeros(number_routes,number_edges);
for m = 1:number_routes
    x_st(m,:) = results{m,6};
end

number_binary_variables = number_routes;
M = 100;
%objective function
z_obj = 1;
y_1_obj = zeros(1,number_edges);
y_2_obj = zeros(1,number_edges);
binary_obj = zeros(1,number_binary_variables);

obj = [z_obj y_1_obj y_2_obj binary_obj];

%first constraint: shortest path longer than any path after upgrading
y_1_first_constraint = red_shsr.*x_st;
y_2_first_constraint = red_hsr.*x_st;
z_first_constraint = ones(number_routes,1);
binary_first_constraint = kron(eye(number_routes),M);
first_constraint_matrix = [z_first_constraint y_1_first_constraint y_2_first_constraint binary_first_constraint];
rhs_first_constraint = t_zero*transpose(x_st);
rhs_first_constraint = transpose(rhs_first_constraint);

%extra binary constraint: b_i1 + ... + b_iK = K-1 for all i
z_binary_constraint = zeros(length(all_K),1);
y_1_binary_constraint = zeros(length(all_K),number_edges);
y_2_binary_constraint = zeros(length(all_K),number_edges);

binary_matrix = zeros(length(all_K), sum(all_K));
row = 1;
column = 1;

for n = 1:length(all_K)
    element = all_K(n);
    counter = 1;
    for operations = 1:element
        if counter <= element
            binary_matrix(row, column) = 1;
            column = column + 1;
            counter = counter + 1;
        end
    end
    row = row + 1;
end
binary_constraint_matrix = [z_binary_constraint y_1_binary_constraint y_2_binary_constraint binary_matrix];
rhs_binary_constraint = transpose(all_K - 1);

%second constraint: edge can only be upgraded once
y_1_second_constraint = eye(number_edges);
y_2_second_constraint = eye(number_edges);
z_second_constraint = zeros(number_edges,1);
binary_second_constraint = zeros(number_edges,number_binary_variables);
second_constraint_matrix = [z_second_constraint y_1_second_constraint y_2_second_constraint binary_second_constraint];
rhs_second_constraint = ones(number_edges,1);

%third constraint: only upgradeable edges are upgraded to shsr
z_third_constraint = zeros(number_edges,1);
y_1_third_constraint = eye(number_edges);
y_2_third_constraint = zeros(number_edges);
binary_third_constraint = zeros(number_edges,number_binary_variables);
third_constraint_matrix = [z_third_constraint y_1_third_constraint y_2_third_constraint binary_third_constraint];
rhs_third_constraint = u_shsr;

%fourth constraint: only upgradeable edges are upgraded to hsr
z_fourth_constraint = zeros(number_edges,1);
y_1_fourth_constraint = zeros(number_edges);
y_2_fourth_constraint = eye(number_edges);
binary_fourth_constraint = zeros(number_edges,number_binary_variables);
fourth_constraint_matrix = [z_fourth_constraint y_1_fourth_constraint y_2_fourth_constraint binary_fourth_constraint];
rhs_fourth_constraint = u_hsr;


%fifth constraint: cost within budget
z_fifth_constraint = 0;
y_1_fifth_constraint = cost_shsr;
y_2_fifth_constraint = cost_hsr;
binary_fifth_constraint = zeros(1,number_binary_variables);
fifth_constraint_matrix = [z_fifth_constraint y_1_fifth_constraint y_2_fifth_constraint binary_fifth_constraint];



constraint_matrix = [first_constraint_matrix; binary_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; fifth_constraint_matrix];
sparse_constraint_matrix = sparse(constraint_matrix);
rhs_till_now = [rhs_first_constraint; rhs_binary_constraint; rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint];

rhs_fifth_constraint = budget;
rhs = [rhs_till_now; rhs_fifth_constraint];
model = [];
model.obj = obj;
model.A = sparse_constraint_matrix;
model.rhs = rhs;
model.sense = [repmat('>',number_routes,1); repmat('=',length(all_K),1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<'];
model.modelsense = 'min';
model.vtype = ['C';repmat('B',2*number_edges,1);repmat('B',number_binary_variables,1)];
params.outputflag = 0;

solution = gurobi(model,params);
optimal_value = solution.objval;

upgrades = solution.x;
upgrades(1) = [];
upgrades_shsr = upgrades(1:number_edges);
upgrades_hsr = upgrades(number_edges+1:2*number_edges);

used_budget = cost_shsr*upgrades_shsr + cost_hsr*upgrades_hsr;

upgrades_shsr_copy = upgrades_shsr;
upgrades_hsr_copy = upgrades_hsr;
nodes = size(hours,1);
upgrade_matrix_hsr = zeros(nodes);
row_counter = 0;
for i = 1:nodes
   row_counter = row_counter + 1;
   for j = 1:nodes
      if j + row_counter > nodes
          break
      end
      upgrade_matrix_hsr(i,j+row_counter) = upgrades_hsr_copy(1);
      upgrades_hsr_copy(1) = [];
   end
end

upgrade_matrix_shsr = zeros(nodes);
row_counter = 0;
for i = 1:nodes
   row_counter = row_counter + 1;
   for j = 1:nodes
      if j + row_counter > nodes
          break
      end
      upgrade_matrix_shsr(i,j+row_counter) = upgrades_shsr_copy(1);
      upgrades_shsr_copy(1) = [];
   end
end

all_routes_selected = zeros(length(all_K),1);
counter = 1;
routes_selected = solution.x(2*number_edges+2:end);
for z = 1:length(all_K)
    routes = all_K(z);
    routes_selected_subset = routes_selected(1:routes);
    [row,~] = find(routes_selected_subset == 0);
    routes_selected(1:routes) = [];
    all_routes_selected(counter,1) = row;
    counter = counter + 1;
end