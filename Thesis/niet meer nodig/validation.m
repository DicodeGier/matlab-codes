addpath("C:\gurobi1001\win64\matlab")

clear
clc


hours = [
  0 2 2 4 ;
  2 0 5 1 ;
  2 5 0 0 ;
  4 1 0 0 ;
];

nodes = 4;
possibilities = nodes*(nodes - 1)/2;
result = cell(possibilities,6);

counter = 1;
for i = 1:nodes
    for j = i+1:nodes
        [cost, route] = dijkstra(hours,i,j);
        result{counter,1} = i;
        result{counter,2} = j;
        result{counter,3} = route;
        result{counter,4} = cost;
        matrix = zeros(nodes);
        for k = 1:length(route)-1
           element_1 = route(k);
           element_2 = route(k+1);
           row_node = min(element_1,element_2);
           column_node = max(element_1,element_2);
           matrix(row_node,column_node) = 1;
        end
        result{counter,5} = matrix;
        indices = triu(ones(nodes),1);
        matrix(~indices) = Inf;
        matrix = matrix';
        matrix = matrix(:);
        matrix(matrix==Inf) = [];
        matrix = matrix';
        result{counter,6} = matrix;
        counter = counter + 1;
        
    end
end

number_routes = size(result,1);
number_edges = length(result{1,6});
x_st = zeros(number_routes,number_edges);
for k = 1:number_routes
    x_st(k,:) = result{k,6};
end

t = [2 2 4 3 5 1 0 3 0 0];
shsr = [2 2 3 2 4 1 0 2 0 0];
red_shsr = t - shsr;
hsr = [1 1 2 1 3 1 0 1 0 0];
red_hsr = t - hsr;
u_shsr = [0;0;1;1;1;0;0;1;0;0];
cost_shsr = 15.*u_shsr;
u_hsr = [1;1;1;1;1;0;0;1;0;0];
cost_hsr = 20.*u_hsr;
budget = 15;

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
fifth_constraint_matrix = [z_fifth_constraint transpose(y_1_fifth_constraint) transpose(y_2_fifth_constraint)];
rhs_fifth_constraint = budget;

constraint_matrix = [first_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; fifth_constraint_matrix];
rhs = [transpose(rhs_first_constraint); rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint; rhs_fifth_constraint];

%gurobi model
model = [];
model.obj = obj;
model.A = sparse(constraint_matrix);
model.rhs = rhs;
model.sense = [repmat('>',number_routes,1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<'];
model.modelsense = 'min';
model.vtype = ['C';repmat('B',2*number_edges,1)];

solution = gurobi(model)

%create matrix with edges that are upgraded for better visualization
upgrades = solution.x;
upgrades(1) = [];
upgrades_shsr = upgrades(1:number_edges);
upgrades_hsr = upgrades(number_edges+1:end);

upgrade_matrix_hsr = zeros(nodes);
row_counter = 0;
for i = 1:nodes
   row_counter = row_counter + 1;
   for j = 1:nodes
      if j + row_counter > nodes
          break
      end
      upgrade_matrix_hsr(i,j+row_counter) = upgrades_hsr(1);
      upgrades_hsr(1) = [];
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
      upgrade_matrix_shsr(i,j+row_counter) = upgrades_shsr(1);
      upgrades_shsr(1) = [];
   end
end

fprintf('maximum shortest path after upgrading is %g hours\n', solution.objval)

upgrade_matrix_shsr
upgrade_matrix_hsr