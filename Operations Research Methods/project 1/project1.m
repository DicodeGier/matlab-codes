clear
clc

%addpath("C:\gurobi1000\win64\matlab")
%nchoosek generates a lot of warnings, surpress them with the following
%line
warning('off')
 
C = load('TSP5.mat');
%use the command below when C is a struct
C = C.C;


%%%parameters for function
k2 = 10000;
k3 = 15000;
E = 20000;
w = ones(size(C,1),1);
%w(10:15,1) = 0

%{
%code for 4d plot
step_x = 10;
step_y = 50;
step_z = 200;
x_values = [];
y_values = [];
z_values = [];
obj_vals = [];
for x_value = 1:step_x:k2
    for y_value = 1:step_y:k3
        for z_value = 1:step_z:E
            x_values(end+1) = x_value;
            y_values(end+1) = y_value;
            z_values(end+1) = z_value;
            [lambda,z_ld,z_lr,check] = tspLD(C,x_value,y_value,z_value,w);
            obj_vals(end+1) = z_ld;
        end
    end
end
scatter3(x_values,y_values,z_values,40,obj_vals,'filled')
xlabel('k2')
ylabel('k3')
zlabel('E')
title("value of objective value under different values for k2, k3 and E")
cb = colorbar;
cb.Label.String = 'value of objective function';
%}


row = 0

%Double for loop for plotting
% i goes over the x-axis
% j goes over the y-axis
x_ax = k2;
y_ax = k3;
% steps = amount of steps on the x a and y axis
% NOTE x_ax and y_ax must be a multiple of steps otherwise it WONT work
steps = 5;
for i = 1:x_ax/steps:x_ax
    column = 0;
    row = row + 1;
    for j = 1:y_ax/steps:y_ax
        column = column + 1;
        %CHANGE THIS LINE WHEN GENERATING A NEW PLOT
        %(which input arguments are i and j
        %ALSO CHANGE THE TITEL AND LABELS FOR THE PLOT BELOW
        [lambda,z_ld,z_lr,check] = tspLD(C,i,j,E,w);
        z(row,column) = z_ld
    end
end
[x,y] = meshgrid(1:x_ax/steps:x_ax,1:y_ax/steps:y_ax);
surf(x,y,z)
title("plot of objective value of Lagrangian dual under different values of k3 and E")
xlabel('k2')
ylabel('k3')
zlabel('objective vale LD')
%}

%[lambda,z_ld,z_lr,check] = tspLD(C,k2,k3,E,w)

function [lambda, z_ld, z_lr,check] = tspLD(C, k2, k3, E, w)
%please set check_bounds = 0 if you do not want to check lower and upper bounds
check_bounds = 0;
%note: output x is of the form [x11 x21 x31 ... ]
%i.e. we strech out the X matrix into a vector and do this by stacking the
%columns of the X matrix
n = size(C,1);

%Enforce S ⊆ W
loop_vec = 1:n;
loop_vec = loop_vec(w == 1);
n_loop_vec = size(loop_vec,2);

%Calculates maximum amount of constraints possible in D (dualized
%polyhedron)
possible_sizes = 4:n_loop_vec-1;
max_E = 0;
for i = possible_sizes
    max_E = max_E + nchoosek(n_loop_vec, i);
end

%Checks wheter user inputted k2, k3 and E are eligible and not larger than
%the possible amount of constraints
if k2 > nchoosek(n,2) || k3 > nchoosek(n_loop_vec ,3) || E > max_E
    lambda = 0;
    z_ld = 0;
    z_lr = 0;
    check = "cannot perform operation with chosen parameters";
else
%start creating the constraints
ones_vector = ones(1,n);
first_constraint = kron(eye(n),ones_vector);
second_constraint = repmat(eye(n),1,n);
constraint_matrix = [first_constraint;second_constraint];

%preallocate rhs, already include rhs of k2 and k3 constraints
rhs = [ones(2*n,1);ones(k2,1);repmat(2,k3,1)];

%generate k2 and k3 constraints
constraint_matrix = constraint_generator(constraint_matrix,n,loop_vec,k2,2);
constraint_matrix = constraint_generator(constraint_matrix,n,loop_vec,k3,3);

%preallocate for dualization
dualized_indices = zeros(E,size(loop_vec,2));
dualized_constraint_matrix = zeros(E,n^2);
dualized_rhs = zeros(E,1);

%For all dualized constraints we try to take all possible constraints for the
%largest subsets of vertices in W 
sizes = zeros(1,E);
iter = 1;
start_point = 1;
while size(sizes,2) <= E
    amount_to_add = min(E - start_point + 1, nchoosek(n_loop_vec,n_loop_vec-iter));
    sizes(1,start_point:start_point+amount_to_add) = n_loop_vec - iter;
    start_point = start_point + amount_to_add;
    iter = iter + 1;
end
sizes = sizes(1,1:E);

%now create constraints to be dualized
for i = 1:E
    if i == 1 | sizes(i-1) ~= sizes(i)
        subhandle = loopchoose(loop_vec,sizes(1,i));
    end
    indices_4 = subhandle();
    dualized_indices(i,1:size(indices_4,2)) = indices_4;
    output_matrix_4 = zeros(n);
    output_matrix_4(indices_4,indices_4) = 1;
    constraint_to_add_4 = output_matrix_4(:);
    dualized_constraint_matrix(i,:) = transpose(constraint_to_add_4);
    dualized_rhs(i,1) = sizes(i)-1;
end

%Presets subgradient algorithm
%NOTE we take a negative lambda
lambda = -1 * ones(E , 1); 
stop = 0;
k = 1;
mu_zero = 1;
alpha = 0.95;
epsilon = 0.01;
threshold = 0.001;
lambdas(:,1) = [lambda];
z_lr = [0];
 
%We can define the lagrangian relaxation outside of the loop with the
%objective function as the exeption as this is the only component changing
%with lambda
main_model.A = sparse(constraint_matrix);
main_model.rhs = rhs;
main_model.sense = [repmat('=',2*n,1);repmat('<',k2+k3,1)];
main_model.modelsense = 'min';
main_model.vtype = repmat('B',size(constraint_matrix,2),1);
params.outputflag = 0;
while stop == 0
    %use gurobi to calculate results
    main_model.obj = transpose(C(:)) - transpose(lambda) * dualized_constraint_matrix;
    result = gurobi(main_model, params);
    
    %add result of current loop iteration
    z_lr(1,end+1) = result.objval+lambda'*dualized_rhs;
    
    %stopping criterion for objective value
    if abs(z_lr(end) - z_lr(end-1)) < threshold
        break
    end
    
    %calculate gamma and determine if algorithm should continue
    gamma_x = dualized_rhs - dualized_constraint_matrix * result.x;
    
    if gamma_x == zeros(E,1)%1 element van lambda of alle elements??
        fprintf('iteration: %g\n',k)
        fprintf('gamma is the zero vector in this iteration: stop algorithm\n')
        break
    end
    
    %determine step size mu
    mu = mu_zero * (alpha^k);
    
    %output for progress loop
    fprintf('iteration: %g\n',k)
    fprintf('stepsize: %g\n',mu)
    fprintf('current solution of Lagrangian relaxation: %g\n',z_lr(1,end))
    fprintf('\n')
    
    %determine lambda_k+1 and determine if algorithm should continue
    lambda_next = min(lambda + mu*gamma_x,0);
    if norm(lambda_next - lambda) > epsilon
        k = k+1;
        lambda = lambda_next;
        lambdas(:,end+1) = lambda;
    else
        break
    end

end

%algorithm has found an optimal solution, now calculate a lower and upper
%bound if desired by user
if check_bounds == 1
params.outputflag = 0;
model_integer = [];
model_integer.obj = C(:);
model_integer.A = sparse([constraint_matrix;dualized_constraint_matrix]);
model_integer.rhs = [rhs;dualized_rhs];
model_integer.sense = [repmat('=',2*n,1);repmat('<',k2+k3+E,1)];
model_integer.modelsense = 'min';
model_integer.vtype = repmat('B',size(constraint_matrix,2),1);
result_integer = gurobi(model_integer,params);

params.outputflag = 0;
model_con = [];
model_con.obj = C(:);
model_con.A = sparse([constraint_matrix;dualized_constraint_matrix]);
model_con.rhs = [rhs;dualized_rhs];
model_con.sense = [repmat('=',2*n,1);repmat('<',k2+k3+E,1)];
model_con.modelsense = 'min';
model_con.lb = zeros(size(constraint_matrix,2),1);
model_con.ub = ones(size(constraint_matrix,2),1);
model_con.vtype = repmat('C',size(constraint_matrix,2),1);
result_con = gurobi(model_con,params);

z_ld = z_lr(end);

lambda
fprintf('optimal value Lagrange relaxation: %g\n',z_ld)
z_lr

check = [result_con.objval z_ld result_integer.objval];
fprintf('lastly, we check whether the found value lies between the lower bound z_lp and upper bound z_ip: \n')
fprintf('lower bound: %g\n',result_con.objval)
fprintf('optimal solution: %g\n',z_ld)
fprintf('upper bound: %g',result_integer.objval)

%if check is not desired
else
z_ld = z_lr(end);

lambda
fprintf('optimal value Lagrange relaxation: %g\n',z_ld)
z_lr

check = [];
check 
end


end
end

function constraint_matrix = constraint_generator(constraint_matrix,n,loop_vec,k,size_set)
    preallocated = zeros(k, n*n);
    if size_set == 2 %here we may take S in V
        subhandle = loopchoose(1:n,size_set);
    else %and here we have to take S in W
        subhandle = loopchoose(loop_vec,size_set);
    end
    to_add = zeros(k,size_set);
    for i = 1:k
        to_add(i,:) = subhandle();
        indices = to_add(i,:);
        output_matrix = zeros(n);
        output_matrix(indices,indices) = 1;
        constraint_to_add = output_matrix(:);
        preallocated(i, :) = transpose(constraint_to_add);
    end
    constraint_matrix = [constraint_matrix; preallocated];

end


%to do for later
%OPSLAAN IN FUNCTIE FILE
%conclusie report
%w veranderen ??
%final output: hoe geef je vectors netter weer
 

%in de gaten houden

