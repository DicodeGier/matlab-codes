clear
clc
addpath("C:\gurobi1001\win64\matlab")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
cost_km_shsr = 16.5;
cost_km_hsr = 25;
budget = 10000;
K = 2;
budgets = 1000:1000:140000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if i == 5 || i == 6 || i == 7 || i == 9 || i == 16 || i == 17 || i == 18 || i == 23 || i == 24 || i == 29 || i == 30 || i == 31 || i == 32 || i == 44
        continue
    end
    for j = i+1:size(hours,2)
        if j == 5 || j == 6 || j == 7 || j == 9 || j == 16 || j == 17 || j == 18 || j == 23 || j == 24 || j == 29 || j == 30 || j == 31 || j == 32 || j == 44
            continue
        end
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
rhs_first_constraint = t_zero*transpose(x_st);
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
rhs_mean_final_constraint = t_zero*transpose(sum(x_st));

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
 