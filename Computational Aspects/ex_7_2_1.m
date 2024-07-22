clear
clc
m = 6;
n = 10;
A = randn(m,n);
b = randn(m,1);

rank_check = rank(A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
part_sol = A\b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_partially = A(:,5:10);
first_part_solution = zeros(4,1);
second_part_solution = A_partially\b;

solution = vertcat(first_part_solution, second_part_solution);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = null(A);
%%x = part_sol + V*alpha = 0 (ideally) <=> V*alpha = -part_sol (niet zelf
%%bedacht)
alpha = V\(-part_sol);
least_sol = part_sol + V*alpha;
given_sol = pinv(A)*b;
least_sol, given_sol %%are the same

