%% unconstrained model
clear
clc

addpath("C:\gurobi1000\win64\matlab")


demand = [335 200 140 440 300 200]
h = 0.3
A = 200

constraint_matrix = [
1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 
0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0
0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0
0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1
1 0 0 0 0 0 -10000 0 0 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 -10000 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 -10000 0 0 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 -10000 0 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0 -10000 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0 0 -10000 0 0 0 0 0 0
]

%x = [q1 ... q6 y1 ... y6 x1 ... x6]
model = []
model.obj = [zeros(1,6) repmat(A,1,6) repmat(h,1,6)]
model.A = sparse(constraint_matrix)
model.rhs = [demand zeros(1,6)]'
model.sense = [repmat('=',1,6) repmat('<',1,6)]'
model.modelsense = 'min'
model.vtype = [repmat('I',1,6) repmat('B',1,6) repmat('I',1,6)]

result = gurobi(model)
fprintf('minimum unconstrained value is %g\n',result.objval)

%% constrained model
clear
clc

addpath("C:\gurobi1000\win64\matlab")


demand = [335 200 140 440 300 200]
h = 0.3
A = 200

constraint_matrix = [
1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 
0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0
0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0
0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0 0
0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1 0
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 -1
1 0 0 0 0 0 -600 0 0 0 0 0 0 0 0 0 0 0
0 1 0 0 0 0 0 -600 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 -600 0 0 0 0 0 0 0 0 0
0 0 0 1 0 0 0 0 0 -400 0 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 0 -200 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0 0 -200 0 0 0 0 0 0
]

%x = [q1 ... q6 y1 ... y6 x1 ... x6]
model = []
model.obj = [zeros(1,6) repmat(A,1,6) repmat(h,1,6)]
model.A = sparse(constraint_matrix)
model.rhs = [demand zeros(1,6)]'
model.sense = [repmat('=',1,6) repmat('<',1,6)]'
model.modelsense = 'min'
model.vtype = [repmat('I',1,6) repmat('B',1,6) repmat('I',1,6)]

result = gurobi(model)
fprintf('minimum constrained value is %g\n',result.objval)
