clear
clc
addpath("C:\gurobi1001\win64\matlab")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
cost_km_shsr = 16.5;
cost_km_hsr = 25;
budget = 10;
K = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hours = [Inf Inf 4   3; 
         Inf Inf 2   5;
         4   2   Inf 6;
         3   5   6   Inf];

nodes = 4;
indices = triu(ones(size(hours,1)),1);
indices = indices > 0;      

t = hours;
t(~indices) = 0;
t = t';
t = t(:);
t(t == 0) = [];
t(t==Inf) = 0;
t = t';

u_shsr = [0; 1; 1; 0; 1; 0];
u_hsr = [0; 1; 1; 0; 1; 1];

cost_shsr = [0 5 4 0 4 0];
cost_hsr = [0 7 8 0 5 3];

shsr = [0 3 2 2 3 6];
hsr = [0 1 1 2 2 4];

red_shsr = t - shsr;
red_hsr = t - hsr;

red_shsr(red_shsr == Inf) = 0;
red_hsr(red_hsr == Inf) = 0;

results = {};
all_K = [];
counter = 1 ;   
for i = 1:size(hours,2)
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
binary_obj = zeros(1,2*number_binary_variables);

obj = [z_obj y_1_obj y_2_obj binary_obj];

%%%not necessary for gurobi but used later on in mean path
%first constraint: shortest path longer than any path after upgrading
y_1_first_constraint = red_shsr.*x_st;
y_2_first_constraint = red_hsr.*x_st;
z_first_constraint = ones(number_routes,1);
binary_first_constraint = kron(eye(number_routes),M);
first_constraint_matrix = [z_first_constraint y_1_first_constraint y_2_first_constraint binary_first_constraint];
rhs_first_constraint = t*transpose(x_st);
rhs_first_constraint = transpose(rhs_first_constraint);

%extra binary constraint: b_i1 + ... + b_iK = K-1 for all i
z_binary_constraint = zeros(length(all_K),1);
y_1_binary_constraint = zeros(length(all_K),number_edges);
y_2_binary_constraint = zeros(length(all_K),number_edges);
z_star_binary_constraint = zeros(length(all_K),number_binary_variables);

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
binary_constraint_matrix = [z_binary_constraint y_1_binary_constraint y_2_binary_constraint binary_matrix z_star_binary_constraint];
rhs_binary_constraint = transpose(all_K - 1);
%%%%%

%constraints for mean path
%1: z_star - 100b_u <= 0
y_1_mean_first = zeros(number_routes,number_edges);
y_2_mean_first = zeros(number_routes,number_edges);
z_mean_first = zeros(number_routes,1);
binary_mean_first = kron(eye(number_routes),-1*M);
z_star_mean_first = eye(number_routes);

mean_first_constraint_matrix = [z_mean_first y_1_mean_first y_2_mean_first binary_mean_first z_star_mean_first];
rhs_mean_first_constraint = zeros(number_routes,1);

%2: z_star >= 0
y_1_mean_second = zeros(number_routes,number_edges);
y_2_mean_second = zeros(number_routes,number_edges);
z_mean_second = zeros(number_routes,1);
binary_mean_second = zeros(number_routes,number_binary_variables);
z_star_mean_second = eye(number_routes);

mean_second_constraint_matrix = [z_mean_second y_1_mean_second y_2_mean_second binary_mean_second z_star_mean_second];
rhs_mean_second_constraint = zeros(number_routes,1);

%3: z_star + y1 * ... + y2 * ... <= t_e * x_e
z_star_mean_third = eye(number_routes);
binary_mean_third = zeros(number_routes);
z_mean_third = zeros(number_routes,1);

mean_third_constraint_matrix = [z_mean_third y_1_first_constraint y_2_first_constraint binary_mean_third z_star_mean_third];
rhs_mean_third_constraint = rhs_first_constraint;

%4: z_star + y1 * ... + y2 * ... + b_u*M >= t_e * x_e - M
z_star_mean_fourth = eye(number_routes);
binary_mean_fourth = kron(eye(number_routes),M);
z_mean_fourth = zeros(number_routes,1);

mean_fourth_constraint_matrix = [z_mean_fourth y_1_first_constraint y_2_first_constraint binary_mean_fourth z_star_mean_fourth];
rhs_mean_fourth_constraint = rhs_first_constraint - M;

%5: final constraint to calculate mean path
z_mean_final = 6;
y_1_mean_final = red_shsr.*sum(x_st);
y_2_mean_final = red_hsr.*sum(x_st);
binary_mean_final = zeros(1,number_binary_variables);
z_star_mean_final = ones(1,number_binary_variables);

mean_final_constraint_matrix = [z_mean_final y_1_mean_final y_2_mean_final binary_mean_final z_star_mean_final];
rhs_mean_final_constraint = t*transpose(sum(x_st));

%second constraint: edge can only be upgraded once
y_1_second_constraint = eye(number_edges);
y_2_second_constraint = eye(number_edges);
z_second_constraint = zeros(number_edges,1);
binary_second_constraint = zeros(number_edges,2*number_binary_variables);
second_constraint_matrix = [z_second_constraint y_1_second_constraint y_2_second_constraint binary_second_constraint];
rhs_second_constraint = ones(number_edges,1);

%third constraint: only upgradeable edges are upgraded to shsr
z_third_constraint = zeros(number_edges,1);
y_1_third_constraint = eye(number_edges);
y_2_third_constraint = zeros(number_edges);
binary_third_constraint = zeros(number_edges,2*number_binary_variables);
third_constraint_matrix = [z_third_constraint y_1_third_constraint y_2_third_constraint binary_third_constraint];
rhs_third_constraint = u_shsr;

%fourth constraint: only upgradeable edges are upgraded to hsr
z_fourth_constraint = zeros(number_edges,1);
y_1_fourth_constraint = zeros(number_edges);
y_2_fourth_constraint = eye(number_edges);
binary_fourth_constraint = zeros(number_edges,2*number_binary_variables);
fourth_constraint_matrix = [z_fourth_constraint y_1_fourth_constraint y_2_fourth_constraint binary_fourth_constraint];
rhs_fourth_constraint = u_hsr;


%fifth constraint: cost within budget
z_fifth_constraint = 0;
y_1_fifth_constraint = cost_shsr;
y_2_fifth_constraint = cost_hsr;
binary_fifth_constraint = zeros(1,2*number_binary_variables);
fifth_constraint_matrix = [z_fifth_constraint y_1_fifth_constraint y_2_fifth_constraint binary_fifth_constraint];



constraint_matrix = [mean_first_constraint_matrix; mean_second_constraint_matrix; mean_third_constraint_matrix; mean_fourth_constraint_matrix; mean_final_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; binary_constraint_matrix; fifth_constraint_matrix];
sparse_constraint_matrix = sparse(constraint_matrix);
rhs_till_now = [rhs_mean_first_constraint; rhs_mean_second_constraint; rhs_mean_third_constraint; rhs_mean_fourth_constraint; rhs_mean_final_constraint; rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint; rhs_binary_constraint];

rhs_fifth_constraint = budget;
rhs = [rhs_till_now; rhs_fifth_constraint];
model = [];
model.obj = obj;
model.A = sparse_constraint_matrix;
model.rhs = rhs;
model.sense = [repmat('<',number_routes,1); repmat('>',number_routes,1); repmat('<',number_routes,1); repmat('>',number_routes,1); '='; repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); repmat('=', size(binary_constraint_matrix,1), 1); '<'];
model.modelsense = 'min';
model.vtype = ['C';repmat('B',2*number_edges,1);repmat('B',number_binary_variables,1);repmat('C',number_binary_variables,1)];
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