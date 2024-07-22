clear
clc

C = [150 175 200 185; 175 200 225 250; 200 225 175 190; 225 250 100000 200];
n = size(C,1);
C = C(:);

first_constraint = repmat(eye(n),1,4);
second_constraint = kron(eye(4),ones(1,4));
third_constraint = zeros(n);
third_constraint(4,3) = 1;
third_constraint = third_constraint(:);
third_constraint = third_constraint';
fourth_constraint = [
     1     0     0     0     0     0     -1     0     0     0     0     0     0     0     0     0
     0     0     -1     0     1     0     0     0     0     0     -1     0     0     0     0     0
     0     0     0     0     0     0     -1     0     1     0     0     0     0     0     -1     0
     0     0     0     0     0     0     0     0     0     0     -1     0     1     0     0     0
     0     0     1     0     -1     0     0     0     0     0     0     0     0     0     0     0
     -1     0     0     0     0     0     1     0     -1     0     0     0     0     0     0     0
     0     0     0     0     -1     0     0     0     0     0     1     0     -1     0     0     0
     0     0     0     0     0     0     0     0     -1     0     0     0     0     0     1     0];

constraint_matrix = [first_constraint;second_constraint;third_constraint;fourth_constraint];

rhs = [ones(2*n,1);0;zeros(8,1)];
%note that x = [x(11) x(21) x(31) x(41) x(12) x(22) ...]
model = [];
model.obj = C';
model.A = sparse(constraint_matrix);
model.rhs = rhs;
model.sense = [repmat('=',2*n,1);'=';repmat('<',8,1)];
model.modelsense = 'min';
model.vtype = repmat('B',n^2,1);
params.outputflag = 0;

result = gurobi(model,params);
x = reshape(result.x,4,4)
obj_val = result.objval

 