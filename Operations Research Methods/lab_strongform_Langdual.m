
%%
%%%%%%%%%%%%
%exercise 1%
%%%%%%%%%%%%
clear
clc

%1
c = [4; 3; 4; 4; 7]
D = [12 13 6 0 1;8 4 9 1 2;2 6 6 0 1;3 5 2 1 8;8 0 5 10 8;2 0 3 4 1]
D_stretched = D(:)

%note that vector of results will be: [y1 ... y5 x11 x21 x31 x41 x51 x12 ..
%x52 x13 .. x53 x14 ... x54 x15 ... x55]

model = []
model.obj = [transpose(c) transpose(D_stretched)]

zeros_block = zeros(size(D,1),size(c,1))
identity_matrices = repmat(eye(size(D,1)),1,size(c,1))
first_constraint = [zeros_block identity_matrices]

minus_ones_vector = -1*ones(size(D,1),1)
first_half = kron(eye(size(c,1)),minus_ones_vector)
second_half = eye(size(D,1)*size(D,2))
second_constraint = [first_half second_half]
all_constraints = [first_constraint; second_constraint]

model.A = sparse(all_constraints)

first_half_rhs = ones(size(D,1),1)
second_half_rhs = zeros(size(D,1)*size(D,2),1)
rhs = [first_half_rhs; second_half_rhs]

model.rhs = rhs
model.sense = [repmat('=',size(D,1),1);repmat('<',size(D,1)*size(D,2),1)]
model.modelsense = 'min'
model.vtype = repmat('B',size(all_constraints,2),1)
solution = gurobi(model)

y = transpose(solution.x(1:5))
final_result = reshape(solution.x(6:end),[6 5])

%%
%2
%first formulation
clear
clc

c = [4; 3; 4; 4; 7]
D = [12 13 6 0 1;8 4 9 1 2;2 6 6 0 1;3 5 2 1 8;8 0 5 10 8;2 0 3 4 1]
D_stretched = D(:)

%note that vector of results will be: [y1 ... y5 x11 x21 x31 x41 x51 x12 ..
%x52 x13 .. x53 x14 ... x54 x15 ... x55]

model = []
model.obj = [transpose(c) transpose(D_stretched)]

zeros_block = zeros(size(D,1),size(c,1))
identity_matrices = repmat(eye(size(D,1)),1,size(c,1))
first_constraint = [zeros_block identity_matrices]

minus_ones_vector = -1*ones(size(D,1),1)
first_half = kron(eye(size(c,1)),minus_ones_vector)
second_half = eye(size(D,1)*size(D,2))
second_constraint = [first_half second_half]
all_constraints = [first_constraint; second_constraint]

model.A = sparse(all_constraints)

first_half_rhs = ones(size(D,1),1)
second_half_rhs = zeros(size(D,1)*size(D,2),1)
rhs = [first_half_rhs; second_half_rhs]

model.rhs = rhs
model.sense = [repmat('=',size(D,1),1);repmat('<',size(D,1)*size(D,2),1)]
model.modelsense = 'min'
model.lb = zeros(size(all_constraints,2),1)
model.ub = ones(size(all_constraints,2),1)
model.vtype = repmat('C',size(all_constraints,2),1)
solution_con = gurobi(model)

y = transpose(solution_con.x(1:5))
final_result = reshape(solution_con.x(6:end),[6 5])
objective_value = solution_con.objval

%LP = ILP because polyhedron is integer

%%
%second formulation
clear
clc

c = [4; 3; 4; 4; 7]
D = [12 13 6 0 1;8 4 9 1 2;2 6 6 0 1;3 5 2 1 8;8 0 5 10 8;2 0 3 4 1]
m = size(D,1)
D_stretched = D(:)

%note that vector of results will be: [y1 ... y5 x11 x21 x31 x41 x51 x12 ..
%x52 x13 .. x53 x14 ... x54 x15 ... x55]

model = []
model.obj = [transpose(c) transpose(D_stretched)]

zeros_block = zeros(size(D,1),size(c,1))
identity_matrices = repmat(eye(size(D,1)),1,size(c,1))
first_constraint = [zeros_block identity_matrices]

minus_identity = -1*eye(size(D,2))
row = repmat(1/m,1,size(D,1))
second_part = kron(eye(size(D,2)),row)
second_constraint = [minus_identity second_part]
all_constraints = [first_constraint; second_constraint]

model.A = sparse(all_constraints)

first_half_rhs = ones(size(D,1),1)
second_half_rhs = zeros(size(D,2),1)
rhs = [first_half_rhs; second_half_rhs]

model.rhs = rhs
model.sense = [repmat('=',size(D,1),1);repmat('<',size(D,2),1)]
model.modelsense = 'min'
model.lb = zeros(size(all_constraints,2),1)
model.ub = ones(size(all_constraints,2),1)
model.vtype = repmat('C',size(all_constraints,2),1)
solution_con = gurobi(model)

y = transpose(solution_con.x(1:5))
final_result = reshape(solution_con.x(6:end),[6 5])
objective_value = solution_con.objval

%%
clear
clc
%%%%%%%%%%%%
%exercise 2%
%%%%%%%%%%%%
model = []
model.obj = [6 5 4]
model.A = sparse([5 4 6;5 5 3])
model.rhs = [10;8]
model.sense = ['<','<']
model.modelsense = 'max'
model.vtype = ['B','B','B']

result = gurobi(model)

%%
%LP relax
clear
clc
model = []
model.obj = [6 5 4]
model.A = sparse([5 4 6;5 5 3])
model.rhs = [10;8]
model.sense = ['<','<']
model.modelsense = 'max'
model.lb = [0,0,0]
model.ub = [1,1,1]
model.vtype = ['C','C','C']

result_relaxed = gurobi(model)

%%
clear
clc

lambda = [1];
k = 1;
epsilon = (10)^(-5);

c = [6 5 4];
A = [5 4 6];
b = [10];
N_constraint = [5 5 3];
rhs_constraint = [8];

alpha = 0.9;
mu_zero = 1;

stop = 0;
while stop == 0
    params.outputflag = 0
    model = [];
    model.obj = c - lambda(end)*A;
    model.A = sparse(N_constraint);
    model.rhs = rhs_constraint;
    model.sense = ['<'];
    model.modelsense = 'max';
    model.vtype = repmat('I',size(N_constraint,2),1);
    result = gurobi(model,params);
    opt_val_lang_dual = result.objval + b*lambda(end);
    x_sol = result.x;
    
    gamma_x = b - A*x_sol;
    if gamma_x == 0
        stop = 1;
    else
        mu = mu_zero * (alpha^k);
        lambda(end+1) = max([lambda(end) - mu*gamma_x,0]);
        if abs(lambda(end) - lambda(end-1)) > epsilon
            k = k + 1;
        else
            stop = 1;
        end
    end
end
lambda_opt = lambda(end)
opt_value = opt_val_lang_dual

%%
clear
clc
%%%%%%%%%%%%
%exercise 3%
%%%%%%%%%%%%
C = [10000 132 217 164 58; 132 10000 290 210 79; 217 290 10000 113 303; 164 210 113 10000 196; 58 79 303 196 10000]
W = [10000 50 65 50 40; 50 10000 30 50 25; 65 30 10000 40 75; 50 50 40 10000 25; 40 25 75 25 10000]
L = 150

%function [LD, lambda] = CMST(C,W,L)
    %note output is of the form [x11 x21 x31 ...]
    n = size(C,1)
    constraint_matrix = [ones(1,n*n)]
    constraint_rhs = [n-1]
    
    for i = 2:n-1
        all_combinations = nchoosek(1:n,i)
        for j = 1:size(all_combinations,1)
            indices = all_combinations(j,:)
            output_matrix = zeros(n)
            output_matrix(indices,indices) = 1
            constraint_to_add = transpose(output_matrix(:))
            constraint_matrix(end+1,:) = constraint_to_add
            constraint_rhs(end+1,1) = i-1
        end
    end
    
    k = 1
    lambda = [1]
    epsilon = 0.05
    mu_zero = 1
    alpha = 0.9
    stop = 0
    z_ld = []
    sol_org = []
    
    while stop == 0
        model = []
        model.obj = C(:) - lambda(end)*-1*W(:)
        model.A = sparse(constraint_matrix)
        model.rhs = constraint_rhs
        model.sense = ['=';repmat('<',size(constraint_matrix,1)-1,1)]
        model.modelsense = 'min'
        model.vtype = repmat('B',size(constraint_matrix,2),1)
        result = gurobi(model)
        
        sol_org(end+1,:) = transpose(C(:))*result.x
        z_ld(end+1,:) = result.objval + lambda(end)*-L
        x_sol = reshape(result.x,[n,n])
        
       gamma_x = -L - transpose(-1*W(:))*result.x;
        if gamma_x == 0
            stop = 1;
        else
            mu = mu_zero * (alpha^k);
            lambda(end+1) = max([lambda(end) + mu*gamma_x,0]);
            if abs(lambda(end) - lambda(end-1)) > epsilon
                k = k + 1;
            else
                stop = 1;
            end
        end
    end 
%end
