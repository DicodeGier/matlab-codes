%%
clear
clc

%%%%%%%%%%%%
%exercise 1%
%%%%%%%%%%%%

%parameters
obj_g = [2 1];
obj_f = [2.2 0.8];
A_g = [1 0];
A_f = [1 0];
b = 24;

feas_g = [3 2; 3 1; 1 2];
rhs_g = [80 50 60];
feas_f = [2 2; 2 5; 2 1];
rhs_f = [70 130 60];

%for the while loop, because reduced cost does not become exactly zero
threshold = 10^(-10);

%initialization, first extreme points
list_vertices = [];
ini_g = intlinprog(-1*obj_g,[1,2],feas_g,rhs_g,[],[],zeros(1,2),[]);
ini_f = intlinprog(-1*obj_f,[1,2],[feas_f;1 0],[rhs_f b-ini_g(1)],[],[],zeros(1,2),[b,inf]);
list_vertices(end+1,:) = transpose(ini_g);
list_vertices(end+1,:) = transpose(ini_f);

helper_matrix = [obj_g;obj_f];
number_lambda = 1;
number_mu = 1;

%set parameter so while loop will start
reduced_cost = 1;

while reduced_cost > 0
    %solve the master problem
    objective = helper_matrix .* list_vertices;
    objective_final = transpose(sum(objective,2));
    constraints = transpose(list_vertices(:,1));
    lambda_constraint = ones(1,number_lambda);
    mu_constraint = ones(1,number_mu);
    lambda_mu_constraint = blkdiag(lambda_constraint,mu_constraint);
    total_constraints = [constraints;lambda_mu_constraint];
    
    model_master = [];
    model_master.obj = objective_final;
    model_master.A = sparse(total_constraints);
    model_master.rhs = [b;1;1];
    model_master.sense = ['<';'=';'='];
    model_master.modelsense = 'max';
    model_master.lb = zeros(1,size(total_constraints,2));
    model_master.vtype = repmat('C',size(total_constraints,2),1);
    
    solution_master = gurobi(model_master);
    
    %solve Giapeto's subproblem
    model_giapeto = [];
    model_giapeto.obj = transpose(obj_g - A_g*solution_master.pi(1));
    model_giapeto.A = sparse(feas_g);
    model_giapeto.rhs = rhs_g;
    model_giapeto.sense = ['<';'<';'<'];
    model_giapeto.modelsense = 'max';
    model_giapeto.lb = zeros(2,1);
    model_giapeto.vtype = repmat('I',2,1);
    
    solution_giapeto = gurobi(model_giapeto);
    reduced_cost_giapeto = solution_giapeto.objval - solution_master.pi(2);
    
    %solve Francesca's subproblem
    model_francesca = [];
    model_francesca.obj = transpose(obj_f - A_f*solution_master.pi(1));
    model_francesca.A = sparse(feas_f);
    model_francesca.rhs = rhs_f;
    model_francesca.sense = ['<';'<';'<'];
    model_francesca.modelsense = 'max';
    model_francesca.lb = zeros(2,1);
    model_francesca.vtype = repmat('I',2,1);
    
    solution_francesca = gurobi(model_francesca);
    reduced_cost_francesca = solution_francesca.objval - solution_master.pi(3);
    
    if reduced_cost_giapeto > threshold | reduced_cost_francesca > threshold
        if reduced_cost_giapeto > reduced_cost_francesca
            list_vertices(end+1,:) = solution_giapeto.x;
            helper_matrix(end+1,:) = obj_g;
            number_lambda = number_lambda + 1;
            reduced_cost = reduced_cost_giapeto;
        else
            list_vertices(end+1,:) = solution_francesca.x;
            helper_matrix(end+1,:) = obj_f;
            number_mu = number_mu + 1;
            reduced_cost = reduced_cost_francesca;
        end
    else
        reduced_cost = 0;
        optimal_value = solution_master.objval
    end    
end

%%
clear
clc
%%%%%%%%%%%%
%exercise 2%
%%%%%%%%%%%%
%parameters
L = [75 85 95];
costs = [90 100 110];
q = [181;139;97;83;41];
s = [17 24 27 32 35];
initial_position = 1; %which roll is used for initial feasible solution
threshold_counter = 3 %if > threshold_counter, while loop will stop

%initialization
initial_matrix = diag(floor(L(initial_position)./s));
minimal_reduced_cost = -10 %set so while loop will start

%create variables for later use
positions = repmat(initial_position,1,size(initial_matrix,1))
cost_vector = repmat(costs(initial_position),1,size(initial_matrix,1))
patterns = []
reduced_costs = []
master_solution_obj_old = 0
counter = 0

while minimal_reduced_cost < 0
    %solve master problem
    master_model = []
    master_model.obj = cost_vector
    master_model.A = sparse(initial_matrix)
    master_model.rhs = q
    master_model.sense = [repmat('>',size(initial_matrix,1),1)]
    master_model.modelsense = 'min'
    master_model.lb = zeros(1,size(initial_matrix,2))
    master_model.vtype = [repmat('C',size(initial_matrix,2),1)]
    master_solution = gurobi(master_model)
    master_solution_obj_new = master_solution.objval
    
    %solve subproblem for every L
    for i = 1:size(L,2)
        sub_model = []
        sub_model.obj = master_solution.pi
        sub_model.A = sparse(s)
        sub_model.rhs = L(i)
        sub_model.sense = ['<']
        sub_model.modelsense = 'max'
        sub_model.lb = zeros(1,size(master_solution.pi,1))
        sub_model.vtype = [repmat('I',size(master_solution.pi,1),1)]
        sub_solution = gurobi(sub_model)
        
        %save results for every subproblem
        patterns(end+1,:) = sub_solution.x
        reduced_costs(end+1,:) = 1 - sub_solution.objval
    end
    %check if while loop is stuck on same objective value
    if master_solution_obj_new == master_solution_obj_old
        counter = counter + 1
        if counter > threshold_counter
            minimal_reduced_cost = 0
        end
    end
    %determine lowest reduced cost with pattern that is not already used
    if minimal_reduced_cost < 0
        master_solution_obj_old = master_solution_obj_new
        alpha = 0
        counter = 0
        while alpha == 0
        minimal_reduced_cost = min(reduced_costs)
        position = find(reduced_costs == minimal_reduced_cost,1)
        column_to_add = transpose(patterns(position,:))
        for k = 1:size(initial_matrix,2)
            check = column_to_add == initial_matrix(:,k)
            summed_up = sum(check)
            if summed_up == size(initial_matrix,1)
                patterns(position,:) = []
                reduced_costs(position,:) = []
                counter = counter + 1
                break
            end
            if k == size(initial_matrix,2)
                alpha = 1
            end
        end
        end
        %now add pattern and save results
        if counter ~= size(L,2)
            positions(:,end+1) = position
            cost_vector(:,end+1) = costs(position)
            initial_matrix(:,end+1) = column_to_add
            patterns = []
            reduced_costs = []
        else
            minimal_reduced_cost = 0
        end
    end
end
solution = master_solution.objval

%%
%question 2-2 --> now with constraints on L
%documentation: see above
L = [75 85 95];
costs = [90 100 110];
q = [181;139;97;83;41];
s = [17 24 27 32 35];
L_constraint = [65;210;54];
initial_position = 2;
threshold_counter = 3

initial_matrix = diag(floor(L(initial_position)./s));
initial_block = zeros(size(L,2),size(initial_matrix,1))
initial_block(initial_position,:) = 1
minimal_reduced_cost = -10

positions = repmat(initial_position,1,size(initial_matrix,1))
cost_vector = repmat(costs(initial_position),1,size(initial_matrix,1))
patterns = []
reduced_costs = []
master_solution_obj_old = 0
counter = 0

while minimal_reduced_cost < 0
    master_model = []
    master_model.obj = cost_vector
    master_model.A = sparse([initial_matrix;initial_block])
    master_model.rhs = [q;L_constraint]
    master_model.sense = [repmat('>',size(s,2),1);repmat('<',size(L,2),1)]
    master_model.modelsense = 'min'
    master_model.lb = zeros(1,size(initial_matrix,2))
    master_model.vtype = [repmat('C',size(initial_matrix,2),1)]
    master_solution = gurobi(master_model)
    master_solution_obj_new = master_solution.objval
    
    for i = 1:size(L,2)
        sub_model = []
        sub_model.obj = master_solution.pi(1:end-size(L,2))
        sub_model.A = sparse(s)
        sub_model.rhs = L(i)
        sub_model.sense = ['<']
        sub_model.modelsense = 'max'
        sub_model.lb = zeros(1,size(master_solution.pi(1:end-size(L,2)),1))
        sub_model.vtype = [repmat('I',size(master_solution.pi(1:end-size(L,2)),1),1)]
        sub_solution = gurobi(sub_model)
        
        patterns(end+1,:) = sub_solution.x
        reduced_costs(end+1,:) = 1 - sub_solution.objval
    end
    if master_solution_obj_new == master_solution_obj_old
        counter = counter + 1
        if counter > threshold_counter
            minimal_reduced_cost = 0
        end
    end
    if minimal_reduced_cost < 0
        master_solution_obj_old = master_solution_obj_new
        alpha = 0
        counter = 0
        while alpha == 0
        minimal_reduced_cost = min(reduced_costs)
        position = find(reduced_costs == minimal_reduced_cost,1)
        column_to_add = transpose(patterns(position,:))
        for k = 1:size(initial_matrix,2)
            check = column_to_add == initial_matrix(:,k)
            summed_up = sum(check)
            if summed_up == size(initial_matrix,1)
                patterns(position,:) = []
                reduced_costs(position,:) = []
                counter = counter + 1
                break
            end
            if k == size(initial_matrix,2)
                alpha = 1
            end
        end
        end
        if counter ~= size(L,2)
            positions(:,end+1) = position
            cost_vector(:,end+1) = costs(position)
            test_vector = zeros(1,size(L,2))
            test_vector(position) = 1
            initial_matrix(:,end+1) = column_to_add
            initial_block(:,end+1) = test_vector
            patterns = []
        reduced_costs = []
        else
            minimal_reduced_cost = 0
        end
    end
end
solution = master_solution.objval

%%
clear
clc
%%%%%%%%%%%%
%exercise 3%
%%%%%%%%%%%%
s = [15, 17, 18, 19, 21, 23, 29, 31, 32, 33, 34, 35]
q = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
L = 87


initial_matrix = eye(size(s,2))
current_solutions = {}
constraint_matrices = {}
final_solutions = {}
x = ini_node_solver(initial_matrix,s,q,L,current_solutions,constraint_matrices,final_solutions)
min_values = []
for j = 1:size(x,2)
    sub_ans = cell2mat(x(1,j))
    min_values(end+1) = sum(sub_ans)
end
final_value = min(min_values)


function x = ini_node_solver(initial_matrix,s,q,L,current_solutions,constraint_matrices,final_solutions)
q1 = q
q2 = q
reduced_cost = -10
    while reduced_cost < 0
        master_model = []
        master_model.obj = ones(1,size(initial_matrix,2))
        master_model.A = sparse(initial_matrix)
        master_model.rhs = q
        master_model.sense = repmat('=',size(initial_matrix,1),1)
        master_model.modelsense = 'min'
        master_model.lb = zeros(size(initial_matrix,2),1)
        master_model.vtype = repmat('C',size(initial_matrix,2),1)
        master_solution = gurobi(master_model)

        sub_problem = []
        sub_problem.obj = transpose(master_solution.pi(1:size(s,2)))
        sub_problem.A = sparse(s)
        sub_problem.rhs = L
        sub_problem.sense = '<'
        sub_problem.modelsense = 'max'
        sub_problem.lb = zeros(size(master_solution.pi(1:size(s,2)),1),1)
        sub_problem.vtype = repmat('B',size(master_solution.pi(1:size(s,2)),1),1)
        sub_solution = gurobi(sub_problem)

        reduced_cost = 1 - sub_solution.objval
        if reduced_cost < 0
            initial_matrix(:,end+1) = [sub_solution.x;zeros(size(initial_matrix,1)-size(sub_solution.x,1),1)]
        end
    end
    x_sol = master_solution.x
    [current_solutions,constraint_matrices,final_solutions] = position_check(initial_matrix,current_solutions,final_solutions, x_sol,constraint_matrices)
    old_min = 0
    while sum(cellfun(@isempty, constraint_matrices)) ~= sum(ones(1,size(constraint_matrices,2)))
        sizes = []
        check = cellfun(@isempty,constraint_matrices)
        for i = 1:size(constraint_matrices,2)
            if check(1,i) == 0
                new_size = size(cell2mat(constraint_matrices(1,i)),1)
                sizes(end+1) = new_size
            end
        end
        new_min = min(sizes)
        if new_min ~= old_min
            old_min = new_min
            q1 = [q1 0]
            q2 = [q2 1]
        end
        [current_solutions,constraint_matrices,final_solutions] = branch_lower(s,q1,L,current_solutions,constraint_matrices,final_solutions)
        [current_solutions,constraint_matrices,final_solutions] = branch_upper(s,q2,L,current_solutions,constraint_matrices,final_solutions)
        positions = find(~cellfun(@isempty,current_solutions))
        current_solutions{positions(1)} = {}
        constraint_matrices{positions(1)} = {}
    end
    x = final_solutions
    
end

function [current_solutions,constraint_matrices,final_solutions] = position_check(initial_matrix,current_solutions,final_solutions, x,constraint_matrices)
    positions = find(x ~= 0 & x ~= 1, 1)
    empty_check = isempty(positions)
    if empty_check == 1
        final_solutions{end+1} = x
    else
        current_solutions{end+1} = x
        position = positions(1)
        row_to_add = zeros(1,size(initial_matrix,2))
        row_to_add(1,position) = 1
        initial_matrix(end+1,:) = row_to_add
        constraint_matrices{end+1} = initial_matrix
    end
end
   
function [current_solutions,constraint_matrices,final_solutions] = branch_lower(s,q1,L,current_solutions,constraint_matrices,final_solutions)
    reduced_cost = -10
    threshold = -10^(-5)
    positions = find(~cellfun(@isempty,constraint_matrices))
    initial_matrix = cell2mat(constraint_matrices(positions(1)))
    while reduced_cost < threshold
        master_model_lower_branch = []
        master_model_lower_branch.obj = ones(1,size(initial_matrix,2))
        master_model_lower_branch.A = sparse(initial_matrix)
        master_model_lower_branch.rhs = q1
        master_model_lower_branch.sense = repmat('=',size(initial_matrix,1),1)
        master_model_lower_branch.modelsense = 'min'
        master_model_lower_branch.lb = zeros(size(initial_matrix,2),1)
        master_model_lower_branch.vtype = repmat('C',size(initial_matrix,2),1)
        master_solution_lower_branch = gurobi(master_model_lower_branch)
        reduced_cost = 0
        
        
        if strcmp('INFEASIBLE',master_solution_lower_branch.status) ~= 1
            sub_problem_lower_branch = []
            sub_problem_lower_branch.obj = transpose(master_solution_lower_branch.pi(1:size(s,2)))
            sub_problem_lower_branch.A = sparse(s)
            sub_problem_lower_branch.rhs = L
            sub_problem_lower_branch.sense = '<'
            sub_problem_lower_branch.modelsense = 'max'
            sub_problem_lower_branch.lb = zeros(size(master_solution_lower_branch.pi(1:size(s,2)),1),1)
            sub_problem_lower_branch.vtype = repmat('B',size(master_solution_lower_branch.pi(1:size(s,2)),1),1)
            sub_solution_lower_branch = gurobi(sub_problem_lower_branch)

            reduced_cost = 1 - sub_solution_lower_branch.objval
            if reduced_cost < threshold & size(initial_matrix,2) <= 40
                initial_matrix(:,end+1) = [sub_solution_lower_branch.x;zeros(size(initial_matrix,1)-size(sub_solution_lower_branch.x,1),1)]
            else
                reduced_cost = 0
            end
        else
            reduced_cost = 0
        end
      
    end
    if strcmp('INFEASIBLE',master_solution_lower_branch.status) ~= 1
        x_sol_lower = master_solution_lower_branch.x
        [current_solutions,constraint_matrices,final_solutions] = position_check(initial_matrix,current_solutions,final_solutions, x_sol_lower,constraint_matrices)    
    else
        current_solutions{end+1} = {}
        constraint_matrices{end+1} = {}
    end
 end

   
function [current_solutions,constraint_matrices,final_solutions] = branch_upper(s,q2,L,current_solutions,constraint_matrices,final_solutions)
    positions = find(~cellfun(@isempty,constraint_matrices))    
    initial_matrix = cell2mat(constraint_matrices(positions(1)))
    threshold = -10^(-5)
    reduced_cost = -10
    while reduced_cost < threshold
        master_model_upper_branch = []
        master_model_upper_branch.obj = ones(1,size(initial_matrix,2))
        master_model_upper_branch.A = sparse(initial_matrix)
        master_model_upper_branch.rhs = q2
        master_model_upper_branch.sense = repmat('=',size(initial_matrix,1),1)
        master_model_upper_branch.modelsense = 'min'
        master_model_upper_branch.lb = zeros(size(initial_matrix,2),1)
        master_model_upper_branch.vtype = repmat('C',size(initial_matrix,2),1)
        master_solution_upper_branch = gurobi(master_model_upper_branch)
        reduced_cost = 0
        
        
        if strcmp('INFEASIBLE',master_solution_upper_branch.status) ~= 1
            sub_problem_upper_branch = []
            sub_problem_upper_branch.obj = transpose(master_solution_upper_branch.pi(1:size(s,2)))
            sub_problem_upper_branch.A = sparse(s)
            sub_problem_upper_branch.rhs = L
            sub_problem_upper_branch.sense = '<'
            sub_problem_upper_branch.modelsense = 'max'
            sub_problem_upper_branch.lb = zeros(size(master_solution_upper_branch.pi(1:size(s,2)),1),1)
            sub_problem_upper_branch.vtype = repmat('B',size(master_solution_upper_branch.pi(1:size(s,2)),1),1)
            sub_solution_upper_branch = gurobi(sub_problem_upper_branch)

            reduced_cost = 1 - sub_solution_upper_branch.objval

            if reduced_cost < threshold & size(initial_matrix,2) <= 40
                initial_matrix(:,end+1) = [sub_solution_upper_branch.x;zeros(size(initial_matrix,1)-size(sub_solution_upper_branch.x,1),1)]
            else
                reduced_cost = 0
            end
        else
            reduced_cost = 0
        end
        
    end
    if strcmp('INFEASIBLE',master_solution_upper_branch.status) ~= 1
        x_sol_upper = master_solution_upper_branch.x
        [current_solutions,constraint_matrices,final_solutions] = position_check(initial_matrix,current_solutions, final_solutions, x_sol_upper,constraint_matrices)
    else
        current_solutions{end+1} = {}
        constraint_matrices{end+1} = {}
    end
end


