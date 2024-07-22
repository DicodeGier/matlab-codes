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


