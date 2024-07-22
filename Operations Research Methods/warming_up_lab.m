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
n = 5;
indices = 1:n;
C_n = indices + indices';
C_n = C_n - diag(diag(C_n) - n^3);

C_n = C_n(:);
%n_opt_assignment = assignment(C_n)

%%c
model = [];
model.obj = C;
length_C = length(C);
identity = eye(sqrt(length_C));
coleq = repmat(identity,1,sqrt(length_C));
A_in = ones(1,sqrt(length_C));
A_in2 = kron(eye(sqrt(length_C)),A_in);
Aeq = [coleq;A_in2];
Aeq = sparse(Aeq);
model.A = Aeq;
beq = ones(2*sqrt(length_C),1);
model.rhs = beq;
model.sense = [repmat('=',2*sqrt(length_C),1)];
model.modelsense = 'min';
lb = zeros(length_C,1);
ub = ones(length_C,1);
model.lb = lb;
model.ub = ub;
model.vtype = repmat('B',size(Aeq,2),1);
%results = gurobi(model)

%%%%%%%%%%%
%exercise 2
%%%%%%%%%%%
C = [10000,132,217,164,58;132,10000,290,210,79;217,290,10000,113,303;164,210,113,10000,196;58,79,303,196,10000];
min_length = MST(C,1)

function length_MST = MST(W,r)
    number_vertices = size(W,1);
    W_helper = W; %copy W for later purposes
    T = [r]; %vector of vertices that are already added to the tree
    length_MST = 0; 
    W_outgoing = W(T,:); %W_outgoing will contain all feasible edges
    counter = 1;
    while counter < number_vertices
        if counter == 1
           %determine minimum length and add it to total length
           min_length = min(min(W_outgoing));
           length_MST = length_MST + min_length;
           %find index and add to T
           [row,col] = find(W == min_length);
           T(end+1) = setdiff(col, T);
           %update W, W_outgoing and counter
           W(:,col) = [];
           W_outgoing(:,col) = [];
           counter = counter + 1;
        else
        %first update W_outgoing with row of point that was added the
        %latest
        W_outgoing(end+1,:) = W(T(end),:);
        %determine min length, find index and update T
        min_length = min(min(W_outgoing));
        length_MST = length_MST + min_length;
        [row,col] = find(W_helper == min_length);
        T(end+1) = setdiff(col, T);
        %find the columns in W and W_outgoing that have to be deleted
        column_tobe_deleted = W_helper(:,T(end));
        for i = 1:size(W,2)
            if W(:,i) == column_tobe_deleted
                W(:,i) = [];
                W_outgoing(:,i) = [];
                break
            end
        end
        counter = counter + 1;
        end
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
    lb = zeros(length_C,1);
    ub = ones(length_C,1);
    x = intlinprog(C,intcon,[],[],Aeq,beq,lb,ub);
    x = reshape(x,sqrt(length_C),sqrt(length_C));
end
