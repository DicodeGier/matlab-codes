addpath("C:\gurobi1000\win64\matlab")

%%
%%%%%%%%%%%%
%excercise 1
%%%%%%%%%%%%
clear
clc
f = [2 1];
A = [3 2; 3 1; 1 2];
b = [80; 50; 60];

model = [];
model.obj = f;
model.A = sparse(A);
model.rhs = b;
model.sense = [repmat('<',size(A,1),1)];
model.modelsense = 'max';
model.vtype = repmat('C', size(A,2),1);
params.outputflag = 0

results = gurobi(model,params)

%1
max_primal = results.objval
max_dual = transpose(b)*results.pi

%2
b_change_finishing = [1;0;0];
max_amount_finishing = transpose(b_change_finishing)*results.pi
%--> not willing to pay for extra hour of finishing

b_change_carpentry = [0;1;0];
max_amount_carpentry = transpose(b_change_carpentry)*results.pi
%--> willing to pay less than 60 cents for extra hour of carpentry

%3
f = [2 1 1.9 2 1.8 2];
A = [3 2 4 3.5 3 2.5; 3 1 2.5 2.6 2.4 2.7; 1 2 1.8 1.6 1.5 2.5];

reduced_cost = f - transpose(results.pi) * A;

%%
%%%%%%%%%%%
%exercise 2
%%%%%%%%%%%
clear
clc
%1a
tic
L = 60;
q = [181;139;97;83;41];
s = [17 24 27 32 35];

value = cutting_stock(L,q,s)


%b
% elapsed time = 0.690419

%c
%almost optimal with (conjecture) a max difference of 2

%2
L = [50,61,63,65,70,80,100,150,200,300]
initial_matrix = diag(floor(L(10)./s))
minimal_reduced_cost = -10

positions = ones(1,size(initial_matrix,1))
patterns = []
reduced_costs = []

while minimal_reduced_cost < 0
    master_model = []
    master_model.obj = ones(1,size(initial_matrix,2))
    master_model.A = sparse(initial_matrix)
    master_model.rhs = q
    master_model.sense = [repmat('>',size(initial_matrix,1),1)]
    master_model.modelsense = 'min'
    master_model.lb = zeros(1,size(initial_matrix,2))
    master_model.vtype = [repmat('C',size(initial_matrix,2),1)]
    master_solution = gurobi(master_model)
    
    for i = 1:size(L,2)
        sub_model = []
        sub_model.obj = -master_solution.pi
        sub_model.A = sparse(s)
        sub_model.rhs = L(i)
        sub_model.sense = ['<']
        sub_model.modelsense = 'min'
        sub_model.lb = zeros(1,size(master_solution.pi,1))
        sub_model.vtype = [repmat('I',size(master_solution.pi,1),1)]
        sub_solution = gurobi(sub_model)
        
        patterns(end+1,:) = sub_solution.x
        reduced_costs(end+1,:) = sub_solution.objval + 1
    end
    minimal_reduced_cost = min(reduced_costs)
    if minimal_reduced_cost < 0
        position = find(reduced_costs == minimal_reduced_cost,1)
        positions(:,end+1) = position
        initial_matrix(:,end+1) = transpose(patterns(position,:))
        patterns = []
        reduced_costs = []
    end
end
solution = master_solution.objval

%3
L = 300;
new_value = cutting_stock(L,q,s)

%%
%%%%%%%%%%%
%exercise 3
%%%%%%%%%%%
clear
clc

L = 60;
N = 300;
s = [17 24 27 32 35];
q = [181 139 97 83 41];
result = naive_method(L,N,s,q)
result_continuous = naive_method_continuous(L,N,s,q)
difference = result-result_continuous

function naive_optimal_continuous = naive_method_continuous(L,N,s,q)
    %create matrix for first constraint
    zeros_block = zeros(length(s),N);
    identity_matrices = repmat(eye(length(s)),1,N);
    first_half_constraints = [zeros_block identity_matrices];
    
    %create matrix for second constraint
    zeros_block_2 = kron(eye(N),-L);
    size_arrays = kron(eye(N),s);
    second_half_constraints = [zeros_block_2 size_arrays];
    
    %combine the two half into one full constraint matrix
    A = [first_half_constraints;second_half_constraints];
    
    %rhs of model
    first_half_rhs = q;
    second_half_rhs = zeros(1,N);
    b = transpose([first_half_rhs second_half_rhs]);
    
    
    model = [];
    model.obj = [ones(1,N) zeros(1,length(s)*N)];
    model.A = sparse(A);
    model.rhs = b;
    model.lb = zeros(1,length(s)*N+N);
    model.ub = [ones(1,N) inf(1,length(s)*N)];
    model.sense = [repmat('>',length(q),1);repmat('<',N,1)];
    model.modelsense = 'min';
    model.vtype = [repmat('C',N,1);repmat('C',length(s)*N,1)];
    results = gurobi(model);
    naive_optimal_continuous = results.objval;
end

function naive_optimal = naive_method(L,N,s,q)
    %create matrix for first constraint
    zeros_block = zeros(length(s),N);
    identity_matrices = repmat(eye(length(s)),1,N);
    first_half_constraints = [zeros_block identity_matrices];
    
    %create matrix for second constraint
    zeros_block_2 = kron(eye(N),-L);
    size_arrays = kron(eye(N),s);
    second_half_constraints = [zeros_block_2 size_arrays];
    
    %combine the two half into one full constraint matrix
    A = [first_half_constraints;second_half_constraints];
    
    %rhs of model
    first_half_rhs = q;
    second_half_rhs = zeros(1,N);
    b = transpose([first_half_rhs second_half_rhs]);
    
    
    model = [];
    model.obj = [ones(1,N) zeros(1,length(s)*N)];
    model.A = sparse(A);
    model.rhs = b;
    model.lb = zeros(1,length(s)*N+N);
    model.sense = [repmat('>',length(q),1);repmat('<',N,1)];
    model.modelsense = 'min';
    model.vtype = [repmat('B',N,1);repmat('I',length(s)*N,1)];
    results = gurobi(model);
    naive_optimal = results.objval;
end

function optimal_value = cutting_stock(L,q,s)
    A_2 = diag(floor(L./s));
    initial_matrix = A_2;
    
    lower_bounds = [];
    upper_bounds = [];
    
    alpha = -1;
    model_helper = [];
    model_helper_2 = [];
    while alpha < 0
        %solving of master problem
        model_helper.A = sparse(initial_matrix);
        f = ones(1,size(model_helper.A,2));
        model_helper.obj = f;
        model_helper.rhs = q;
        model_helper.sense = [repmat('>',size(model_helper.A,1),1)];
        model_helper.modelsense = 'min';
        model_helper.vytpe = repmat('C', size(model_helper.A,2),1);
        params.outputflag = 0;
        mu = gurobi(model_helper,params);
        
        %solving of subproblem
        model_helper_2.obj = transpose(mu.pi);
        model_helper_2.A = sparse(s);
        model_helper_2.rhs = L;
        model_helper_2.sense = [repmat('<',size(model_helper_2.A,1),1)];
        model_helper_2.modelsense = 'max';
        model_helper_2.vtype = repmat('I', size(model_helper_2.A,2),1);
        params.outputflag = 0
        alpha_struct = gurobi(model_helper_2,params);
        alpha = 1 - alpha_struct.objval;
        
        %updating of bounds
        lower_bounds(end+1,1) = mu.objval/(1-alpha);
        upper_bounds(end+1,1) = mu.objval;
        
        %only update matrix A if alpha is smaller than 0
        if alpha < 0
            initial_matrix(:,end+1) = alpha_struct.x;
        end
    end
    
    %plot bounds
    optimal_value = mu.objval
    plot(lower_bounds)
    hold on
    plot(upper_bounds)
    xlabel('iteration')
    ylabel('value')
    title('bounds per iteration')
    legend('lower bound', 'upper bound')
    hold off
end
