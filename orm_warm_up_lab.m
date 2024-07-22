clear
clc

%%%%%%%%%%%%
%excercise 1
%%%%%%%%%%%%

%%a
C = [10000,132,217,164,58;132,10000,290,210,79;217,290,10000,113,303;164,210,113,10000,196;58,79,303,196,10000];
C = C(:);
%opt_assignment = assignment(C)

%%b
n = 4;
indices = 1:n;
C_n = indices + indices';
C_n = C_n - diag(diag(C_n) - n^3);

C_n = C_n(:);
n_opt_assignment = assignment(C_n)

%%c
model = [];
model.obj = C;
length_C = length(C);
identity = eye(sqrt(length_C));
coleq = repmat(identity,1,sqrt(length_C));
A_in = ones(1,sqrt(length_C));
A_in2 = kron(eye(sqrt(length_C)),A_in);
Aeq = [coleq;A_in2];
model.A = sparse(Aeq);
beq = ones(2*sqrt(length_C),1);
model.rhs = beq;
model.sense = [repmat('=',2*sqrt(length_C),1)];
model.modelsense = 'min';
lb = zeros(length_C);
ub = ones(length_C);
model.lb = lb;
model.ub = ub;
model.vtype = repmat('I',2*sqrt(length_C),1);
results = gurobi(model)


%%%%%%%%%%%
%exercise 2
%%%%%%%%%%%
min_length = MST(C,1)

function length_MST = MST(W,r)
    number_vertices = size(W,1)
    whole_set = 1:number_vertices
    T = [r]
    W_outgoing = W(T,1:end~=T)
    while whole_set ~= sort(T)
        min_length = min(min(W_outgoing))
        length_MST = length_MST + min_length
        [row,col] = find(W == min_length)
        T(end+1) = col

    end
end





function x = assignment(C)
    length_C = length(C);
    intcon = 1:length_C;
    identity = eye(sqrt(length_C));
    coleq = repmat(identity,1,sqrt(length_C));
    A_in = ones(1,sqrt(length_C));
    A_in2 = kron(eye(sqrt(length_C)),A_in);
    Aeq = [coleq;A_in2];
    beq = ones(2*sqrt(length_C),1);
    lb = zeros(length_C);
    ub = ones(length_C);
    x = intlinprog(C,intcon,[],[],Aeq,beq,lb,ub);
    x = reshape(x,sqrt(length_C),sqrt(length_C));
end