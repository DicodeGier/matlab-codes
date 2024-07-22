addpath("C:\gurobi1001\win64\matlab")

clear
clc

hours = [
  0 2 2 4;
  2 0 5 1;
  2 5 0 3;
  4 1 3 0
];

nodes = 4;
possibilities = nodes*(nodes - 1)/2;
result = cell(possibilities,6);

counter = 1;
for i = 1:nodes
    for j = i+1:nodes
        [cost, route] = dijkstra(hours,i,j);
        route = flip(route,2);
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
length_route = length(result{1,6});
x_st = zeros(number_routes,length_route);
for k = 1:number_routes
    x_st(k,:) = result{k,6};
end

t = [2 2 4 5 1 3];
shsr = [2 2 3 4 1 2];
red_shsr = t - shsr;
hsr = [1 1 2 3 1  1];
red_hsr = t - hsr;
u_shsr = [0;0;1;1;0;1];
cost_shsr = 15.*u_shsr;
u_hsr = [1;1;1;1;0;1];
cost_hsr = 20.*u_hsr;
budget = 14;

z_obj = 1;
y_1_obj = zeros(1,length_route);
y_2_obj = zeros(1,length_route);

obj = [z_obj y_1_obj y_2_obj];

y_1_first_constraint = red_shsr.*x_st;
y_2_first_constraint = red_hsr.*x_st;
z_first_constraint = ones(number_routes,1);
first_constraint_matrix = [z_first_constraint y_1_first_constraint y_2_first_constraint];
rhs_first_constraint = t*transpose(x_st);

y_1_second_constraint = eye(number_routes);
y_2_second_constraint = eye(number_routes);
z_second_constraint = zeros(number_routes,1);
second_constraint_matrix = [z_second_constraint y_1_second_constraint y_2_second_constraint];
rhs_second_constraint = ones(number_routes,1);

z_third_constraint = zeros(number_routes,1);
y_1_third_constraint = eye(number_routes);
y_2_third_constraint = zeros(number_routes);
third_constraint_matrix = [z_third_constraint y_1_third_constraint y_2_third_constraint];
rhs_third_constraint = u_shsr;

z_fourth_constraint = zeros(number_routes,1);
y_1_fourth_constraint = zeros(number_routes);
y_2_fourth_constraint = eye(number_routes);
fourth_constraint_matrix = [z_fourth_constraint y_1_fourth_constraint y_2_fourth_constraint];
rhs_fourth_constraint = u_hsr;

z_fifth_constraint = 0;
y_1_fifth_constraint = cost_shsr;
y_2_fifth_constraint = cost_hsr;
fifth_constraint_matrix = [z_fifth_constraint transpose(y_1_fifth_constraint) transpose(y_2_fifth_constraint)];
rhs_fifth_constraint = budget;

constraint_matrix = [first_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; fifth_constraint_matrix];
rhs = [transpose(rhs_first_constraint); rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint; rhs_fifth_constraint];


model = [];
model.obj = obj;
model.A = sparse(constraint_matrix);
model.rhs = rhs;
model.sense = [repmat('>',number_routes,1); repmat('<',length_route,1); repmat('<',length_route,1); repmat('<', length_route,1); '<'];
model.modelsense = 'min';
model.vtype = ['I';repmat('B',2*length_route,1)];

solution = gurobi(model)

upgrades = solution.x;
upgrades(1) = [];
upgrades_shsr = upgrades(1:length_route);
upgrades_hsr = upgrades(length_route+1:end);

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

upgrade_matrix_shsr
upgrade_matrix_hsr

