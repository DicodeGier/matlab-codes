function [lambda, zld, zlr] = tspLD(C, K2, K3, E, w)
%surpress warnings from nchoosek
warning('off')
%please set check_bounds = 0 if you do not want to check lower and upper bounds
%please be advised, check_bounds on TSP5 and TSP6 causes memory overflow
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
if K2 > nchoosek(n,2) || K3 > nchoosek(n_loop_vec ,3) || E > max_E
    lambda = 0;
    zld = 0;
    zlr = 0;
    fprintf("cannot perform operation with chosen parameters")
else
%start creating the constraints
ones_vector = ones(1,n);
first_constraint = kron(eye(n),ones_vector);
second_constraint = repmat(eye(n),1,n);
constraint_matrix = [first_constraint;second_constraint];

%preallocate rhs, already include rhs of k2 and k3 constraints
rhs = [ones(2*n,1);ones(K2,1);repmat(2,K3,1)];

%generate k2 and k3 constraints
constraint_matrix = constraint_generator(constraint_matrix,n,loop_vec,K2,2);
constraint_matrix = constraint_generator(constraint_matrix,n,loop_vec,K3,3);

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
if n<20
    penalty_percentage = 0.15;
elseif n<60
    penalty_percentage = 0;
elseif n < 150
    penalty_percentage = 0.5;
elseif n<250
    penalty_percentage = 0.7;
elseif n < 350
    penalty_percentage = 0.8;
elseif n < 450
    penalty_percentage = 0.9;
else
    penalty_percentage = 1;

end

alpha = 0.99 - penalty_percentage*0.09;
epsilon = 0.001 + penalty_percentage*0.009;
threshold = 0.001 + penalty_percentage*0.2;
lambdas(:,1) = [lambda];
zlr = [0];
 
%We can define the lagrangian relaxation outside of the loop with the
%objective function as the exeption as this is the only component changing
%with lambda
main_model.A = sparse(constraint_matrix);
main_model.rhs = rhs;
main_model.sense = [repmat('=',2*n,1);repmat('<',K2+K3,1)];
main_model.modelsense = 'min';
main_model.vtype = repmat('B',size(constraint_matrix,2),1);
params.outputflag = 0;
while stop == 0
    %use gurobi to calculate results
    main_model.obj = transpose(C(:)) - transpose(lambda) * dualized_constraint_matrix;
    result = gurobi(main_model, params);
    
    %add result of current loop iteration
    zlr(1,end+1) = result.objval+lambda'*dualized_rhs;
    
    %stopping criterion for objective value
    if abs(zlr(end) - zlr(end-1)) < threshold
        fprintf('iteration: %g\n',k)
        fprintf('algorithm stopped because of threshold for objval\n')
        break
    end
    
    %calculate gamma and determine if algorithm should continue
    gamma_x = dualized_rhs - dualized_constraint_matrix * result.x;
    
    if gamma_x == zeros(E,1)
        fprintf('iteration: %g\n',k)
        fprintf('gamma is the zero vector in this iteration: stop algorithm\n')
        break
    end
    
    %determine step size mu
    mu = mu_zero * (alpha^k);
    
    %output for progress loop
    fprintf('iteration: %g\n',k)
    fprintf('stepsize calculated in this iteration: %g\n',mu)
    fprintf('current solution of Lagrangian relaxation: %g\n',zlr(1,end))
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
model_integer.sense = [repmat('=',2*n,1);repmat('<',K2+K3+E,1)];
model_integer.modelsense = 'min';
model_integer.vtype = repmat('B',size(constraint_matrix,2),1);
result_integer = gurobi(model_integer,params);

params.outputflag = 0;
model_con = [];
model_con.obj = C(:);
model_con.A = sparse([constraint_matrix;dualized_constraint_matrix]);
model_con.rhs = [rhs;dualized_rhs];
model_con.sense = [repmat('=',2*n,1);repmat('<',K2+K3+E,1)];
model_con.modelsense = 'min';
model_con.lb = zeros(size(constraint_matrix,2),1);
model_con.ub = ones(size(constraint_matrix,2),1);
model_con.vtype = repmat('C',size(constraint_matrix,2),1);
result_con = gurobi(model_con,params);

zld = zlr(end);

fprintf('optimal value Lagrange relaxation: %g\n',zld)

fprintf('lastly, we check whether the found value lies between the lower bound z_lp and upper bound z_ip: \n')
fprintf('lower bound: %g\n',result_con.objval)
fprintf('optimal solution: %g\n',zld)
fprintf('upper bound: %g',result_integer.objval)

%if check is not desired
else
zld = zlr(end);

fprintf('optimal value Lagrange relaxation: %g\n',zld) 
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


