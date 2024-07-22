clear
clc

%exercise 1
c1 = [2 1]
c2 = [2.2 0.8]
A1 = [1 0]
b = 24
A2 = [1 0]
u = [0 0; 0 30; 8 26; 16.6 0]
v = [0 0; 0 26; 15 20; 25 10; 30 0]

objective_part_1 = sum(repmat(c1,size(u,1),1).*u,2)
objective_part_2 = sum(repmat(c2,size(v,1),1).*v,2)
objective = [objective_part_1;objective_part_2]

A_1 = repmat(A1,size(u,1),1).*u
A_1 = A_1(:,1)
A_2 = repmat(A2,size(v,1),1).*v
A_2 = A_2(:,1)
A = [transpose(A_1) transpose(A_2)]
Aeq = [1 1 1 1 0 0 0 0 0; 0 0 0 0 1 1 1 1 1]
beq = [1 1]

x = intlinprog(-1*objective,[],A,b,Aeq,beq,zeros(size(A,2)),[])
objval = transpose(x)*objective