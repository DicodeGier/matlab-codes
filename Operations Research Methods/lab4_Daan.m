% Daan van Turnhout
addpath("C:\gurobi1000\win64\matlab")
addpath("C:\gurobi1001\win64\matlab")

%% Exercise 1
clear
clc

c = [4,3,4,4,7];
D = [12,13,6,0,1;
     8,4,9,1,2;
     2,6,6,0,1;
     3,5,2,1,8;
     8,0,5,10,8;
     2,0,3,4,1];

% (a) use F1 formulation
[x,objval] = F1(c,D,'I');
% (b) F1 relaxation
[x_F1,objval_F1] = F1(c,D,'C');
% (b) F2 relaxation
[x_F2,objval_F2] = F2(c,D);

%% Exercise 2
clear
clc

c = [6,5,4];
A = [5,4,6; 5,5,3];
b = [10;8];

% (a)
% IP
params.outputflag = 0;
model = [];
model.obj = c;
model.A = sparse(A);
model.rhs = b;
model.modelsense = 'max';
model.lb = zeros(size(A,2),1);
model.ub = ones(size(A,2),1);
model.vtype = repmat('I', size(A,2),1);
model.sense = repmat('<',size(A,1),1); 
result_IP = gurobi(model,params);
objval_IP = result_IP.objval

% LP
model.vtype = repmat('C', size(A,2),1);
result_LP = gurobi(model,params);
objval_LP = result_LP.objval

% (b)
[lambda, objval_subgradient] = subgradient()

%% Exercise 3
clear
clc



%% Functions
function [x,objval] = F1(c,D,vtype)
    params.outputflag = 0;
    
    n = size(D,2);
    m = size(D,1);
    
    % variables: [n*m x variables, m y variables]
    D = D';
    c_1 = D(:);
    c_2 = c(:);
    c = [c_1;c_2];
    % constraint 1; only 1 assignment per colomn
    block(1:n) = 1;
    A_1 = kron(eye(m),block);
    A_1 = [A_1, zeros(m,n)];
    % constraint 2; xij<yj = xij - yj = 0
    block = repmat(-1* eye(n),m,1);
    A_2 = eye(n*m);
    A_2 = [A_2, block];
    A = [A_1;A_2];
    % RHS
    b_1 = ones(size(A_1,1),1);
    b_2 = zeros(size(A_2,1),1);
    b = [b_1;b_2];
    
    model = [];
    model.obj = c;
    model.A = sparse(A);
    model.rhs = b;
    model.modelsense = 'min';
    model.lb = zeros(size(A,2),1);
    model.ub = ones(size(A,2),1);
    model.vtype = repmat(vtype, size(A,2),1);
    model.sense = [repmat('=',size(A_1,1),1); repmat('<',size(A_2,1) ,1)]; 
    result = gurobi(model,params);
    
    x = result.x;
    objval = result.objval;
end

function [x,objval] = F2(c,D)
    params.outputflag = 0;
    
    n = size(D,2);
    m = size(D,1);
    
    % variables: [n*m x variables, m y variables]
    D = D';
    c_1 = D(:);
    c_2 = c(:);
    c = [c_1;c_2];
    % constraint 1; only 1 assignment per colomn
    block(1:n) = 1;
    A_1 = kron(eye(m),block);
    A_1 = [A_1, zeros(m,n)];
    % constraint 2; sum xij<myj = sum xij - myj = 0
    A_2 = repmat(eye(n),1,m);
    A_2 = [A_2, -m* eye(n)];
    A = [A_1;A_2];
    % RHS
    b_1 = ones(size(A_1,1),1);
    b_2 = zeros(size(A_2,1),1);
    b = [b_1;b_2];
    
    model = [];
    model.obj = c;
    model.A = sparse(A);
    model.rhs = b;
    model.modelsense = 'min';
    model.lb = zeros(size(A,2),1);
    model.ub = ones(size(A,2),1);
    model.vtype = repmat('C', size(A,2),1);
    model.sense = [repmat('=',size(A_1,1),1); repmat('<',size(A_2,1) ,1)]; 
    result = gurobi(model,params);
    
    x = result.x;
    objval = result.objval;
end

function [lambda,objval] = subgradient()
    params.outputflag = 0;
    
    c = [6,5,4];
    % dualize first row of A
    D = [5,5,3];
    d = 8;
    A = [5,4,6];
    b = 10;
    
    % init params
    threshold = 10^(-5);
    k = 1;
    lambda = [1];
    stop = 0;
    mu_zero = 1;
    alpha = 0.9; % should be 'sufficiently high' 
 
    while stop == 0
        % objective
        obj = c - lambda(k) .*A;
        
        model = [];
        model.obj = obj;
        model.A = sparse(D);
        model.rhs = d;
        model.modelsense = 'max';
        model.vtype = repmat('I', size(D,2),1);
        model.sense = repmat('<',size(D,1),1); 
        result_LR = gurobi(model,params);
        z_LR = result_LR.objval + b'.*lambda(k);

        gamma = b - A*result_LR.x;
        if gamma == 0 % lambda(k) is opt
            stop = 1;
        else
            mu = mu_zero * (alpha^k);
            lambda(k+1) = max( lambda(k) - mu*gamma  ,0);
            if abs(lambda(k+1) - lambda(k)) > threshold
                k = k + 1;
            else
                stop = 1; %lambda(k+1) is opt
            end
        end
    end
    lambda = lambda(k+1);
    objval = z_LR;
end