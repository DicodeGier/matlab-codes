clear
clc
close all

profile on

rng(1)

% Initialization:
M = load('inst3.mat').M; %[3,6,4;2,1,5];%load('inst3.mat').M;
K_max = 8; % Max # elements on right hand side
W_max = 100; % Max sum of weights
P = 1000; % # initial solutions
l = 100; % size tournament mating pool
M_size = 50; % size mating pool
l2 = 5; % size tournament parent selection
maxit = 1000; % stopping criterion: maximum of 1000 iterations
sameit = 50; % stopping criterion: stop if after sameit iterations global opt has not changed
T = 100; % stopping criterion: time
mutation_prob = 0.2;


[best_cost,global_opt,B,S,k,w] = ARCDP2(M,K_max,W_max,P,M_size,l,l2,maxit,sameit,T,mutation_prob);
profile viewer
function [best_cost,global_opt,B,S,k,w] = ARCDP2(A,K_max,W_max,P,M_size,l,l2,maxit,sameit,T,mutation_prob)
    best_cost = zeros(1,maxit);
    global_opt = Inf;
    it = 1;
    print = 1;
    totaltime = 0;
    iteration_same = 0;
    [initial_solution,weights,objective_values] = initial_solutions(A,K_max,W_max,P);
    minimum = min(objective_values);
    best_cost(1,it) = minimum;
    if minimum < global_opt
        global_opt = minimum;
        idx = find(objective_values == global_opt);
        idx = idx(1);
        B = initial_solution(idx,:);
    end
    fprintf('----------- %g -----------\n',it)
    fprintf('best cost : %g\n',minimum)
    fprintf('global opt: %g\n',global_opt)
    fprintf('--------------------------\n')
    fprintf('\n')
    while it < maxit && global_opt ~= 0 && print == 1 && totaltime < T %print == 1 is stopping criterion for sameit: if it > sameit => print = 0 so the new iteration will not be printed in the command window
        tstart = tic;
        it = it+1;
        [mating_pool,mating_pool_obj_values] = mating_pool_creater(objective_values,M_size,l);
        parents = parents_selection(mating_pool,mating_pool_obj_values,P,l2);
        [children,weights,objective_values] = offspring_creation(A,initial_solution,weights,P,parents,K_max,mutation_prob);
        initial_solution = children;
        minimum = min(objective_values);
        best_cost(1,it) = minimum;
        if minimum < global_opt
            global_opt = minimum;
            idx = find(objective_values == global_opt);
            idx = idx(1);
            B = initial_solution(idx,:);
        else
            iteration_same = iteration_same + 1;
            if iteration_same > sameit
                print = 0;
            end
        end
        tend = toc(tstart);
        totaltime = totaltime + tend;
        if print == 1
            fprintf('----------- %g -----------\n',it)
            fprintf('best cost : %g\n',minimum)
            fprintf('global opt: %g\n',global_opt)
            fprintf('total time: %g\n',totaltime)
            fprintf('--------------------------\n')
            fprintf('\n')
        end
    end
    fprintf('global optimum is %g\n',global_opt)
    B = B(1,~cellfun(@isempty,B));
    w = cell2mat(cellfun(@(x) max(max(x)),B,'UniformOutput',false));
    B = cellfun(@(x) x./max(max(max(x)),1),B,'UniformOutput',false);
    k = size(B,2);
    total = zeros(size(A,1),size(A,2));
    for i = 1:k
        matrix = cell2mat(B(1,i));
        total = total + matrix.*w(i);
    end
    S = total - A;
    plot(1:it,best_cost(1:it))
end

function [children,new_weights,objective_values] = offspring_creation(A,initial_solution,weights,P,parents,K_max,mutation_prob)
    num_matrices_vector = sum(~cellfun(@isempty,initial_solution),2); % Find # of matrices of each solution
    rows = size(initial_solution{1,1},1);
    columns = size(initial_solution{1,1},2);
    objective_values = zeros(P,1);
    children = cell(P,K_max);
    new_weights = zeros(P,K_max);
    for row = 1:round(P/2) 
        couple = parents(row,:);
        couple = sort(couple);
        parent1 = couple(1);
        parent2 = couple(2);
        num_matrices = num_matrices_vector(couple,:); % Find # matrices of both parents

        if num_matrices(1) > num_matrices(2)
            [parent2,parent1] = deal(parent1,parent2);
        end
        
        num_matrices = num_matrices_vector([parent1,parent2],:);
        
        parent_cell1 = initial_solution(parent1,:);
        parent_cell2 = initial_solution(parent2,:);

        weight_parent1 = weights(parent1,:);
        weight_parent2 = weights(parent2,:);
 
        max_matrices = max(num_matrices);
        min_matrices = min(num_matrices);

        % Creating offspring:
        j = 1;
        for i = 1:max_matrices
            matrix1 = cell2mat(parent_cell1(1,j));
            matrix2 = cell2mat(parent_cell2(1,i));
            offspring_matrix1 = matrix1;
            offspring_matrix2 = matrix2;

            % determine the weight of both matrices
            weight1 = weight_parent1(j);
            weight2 = weight_parent2(i);
            childweight = randi([min(weight1,weight2),max(weight1,weight2)]);%(weight1 + weight2)/2;
            if mod(childweight,1) ~= 0
                childweight = round(childweight) - 1;
            end
            
            u = rand;
            if u <= mutation_prob
                if childweight == 1
                    childweight = childweight + 1;
                else
                    childweight = childweight - 1;
                end
            end

            % Swap a random row of the matrices to create to new offspring matrices
            
            row_tbs = randperm(rows,round(rows/2));
            offspring_matrix1(row_tbs,:) = matrix2(row_tbs,:);
            offspring_matrix2(row_tbs,:) = matrix1(row_tbs,:);
            
            
            % Add found matrices to the child 3D matrices
            children{2*row-1,i} = offspring_matrix1;
            new_weights(2*row-1,i) = childweight;
            children{2*row,i} = offspring_matrix2;
            new_weights(2*row,i) = childweight;

            if j < min_matrices % Makes sure that if we reach the end of parent 1's
                                % matrices we go back to the first matrix
                j = j+1;
            else
                j = 1;
            end
        end
        child1_3d = reshape(cell2mat(children(2*row-1,:)),rows,[],max_matrices);
        weight = new_weights(2*row-1,:);
        weight = weight(weight > 0)';
        child1_3d_weighted = child1_3d.*reshape(weight,1,1,[]);
        objective_values(2*row-1,1) = sum(sum((sum(child1_3d_weighted,3) - A).^2));
        child2_3d = reshape(cell2mat(children(2*row,:)),rows,[],max_matrices);
        weight = new_weights(2*row,:);
        weight = weight(weight > 0)';
        child2_3d_weighted = child2_3d.*reshape(weight,1,1,[]);
        objective_values(2*row,1) = sum(sum((sum(child2_3d_weighted,3) - A).^2));
    end
end

function parents = parents_selection(mating_pool,mating_pool_obj_values,P,l2)
    % Selects parent couples from the mating pool and stores these pairs in
    % a matrix. Takes as input variables a row vector as mating pool, # of 
    % initial solutions and size of the tournament.
    M_size = length(mating_pool);
    parents = zeros(round(P/2),2);
    
    i = 1;
    while i <= round(P/2)
        indices = randperm(M_size,l2); % Choose l2 distinct random indices between 1 and M_size
        values = mating_pool_obj_values(indices); %Find corresponding values of selected indices
        [~,row] = min(values);
        mating_pool_obj_values(indices(row)) = Inf; %Select index with the minimal value and set its value to +Inf
        indices2 = randperm(M_size,l2);
        values2 = mating_pool_obj_values(indices2);
        [~,row2] = min(values2);
        % Append the indices of both parents to the matrix "parents"
        parents(i,1) = mating_pool(indices(row));
        parents(i,2) = mating_pool(indices2(row2));
        i = i+1;
    end
 end



 function [mating_pool,mating_pool_obj_values] = mating_pool_creater(objective_values,M_size,l)
    % Creates a mating pool of size M by tournament selection. Takes as
    % input variables the objective values of the solutions,the size of the
    % mating pool and the size of the tournament.
    P = length(objective_values);
    mating_pool = zeros(1,M_size);
    
    i = 1;
    mating_pool_obj_values = zeros(1,M_size);
    while i <= M_size
        indices = randperm(P,l); % Select l district random indices between 1 and P
        values = objective_values(indices); %Find corresponding objective values of selected indices
        [mini,row] = min(values);
        mating_pool_obj_values(1,i) = mini;
        mating_pool(1,i) = indices(row); % Add the index of the solution in the mating pool with the minimal value
        i = i+1;
    end
 end


 function [all_solutions,weights,objective_values] = initial_solutions(A,K_max,W_max,P)
    % Generates initial solutions to the problem. Takes as input variables
    % the left-hand side matrix A, max number of matrices, max sum of weights
    % and the desired number of intial solutions.
    
    % Generate a random vector that determines the number of matrices on
    % the right-hand side for each individual solutions and sort this
    % vector in appending order:
    number_matrices = randi([1,K_max],P,1);
    number_matrices = sort(number_matrices);
        
    rows = size(A,1);
    columns = size(A,2);
    all_solutions = cell(P,K_max);
    weights = zeros(P,K_max);
    objective_values = zeros(P,1);
    
    for j = 1:P
        summed_up = zeros(rows,columns);
        elements = number_matrices(j);
        weight_left = W_max;
        max_weight = W_max-elements+1;
        matrices_left = elements;
        % Loop for creating matrices & corresponding weights:
        for k = 1:elements
            matrix = zeros(rows,columns);
            % Loop for creating a row in a matrix:
            alpha = 1;
            while alpha == 1
                for i = 1:rows
                    integer1 = randi(columns);
                    integer2 = randi(columns);
                    first_integer = min(integer1,integer2);
                    second_integer = max(integer1,integer2);
                    first_array = zeros(1,first_integer);
                    second_array = ones(1,second_integer-first_integer);
                    third_array = zeros(1,columns - second_integer);
                    row = [first_array second_array third_array];
                    matrix(i,:) = row;
                end
                if all(matrix == 0,'all') ~= 1
                    alpha = 0;
                end     
            end
            % Selecting a weight for the created matrix:
            weight = randi(max_weight);
            weights(j,k) = weight;
            max_weight = weight_left - weight - matrices_left + 1 + 1;
            weight_left = weight_left - weight;
            matrices_left = matrices_left - 1;
            
            %matrix = matrix.*weight; % Multiply matrix by its weight to easily store the weight
            summed_up = summed_up + weight.*matrix; % Add new found weight*matrix to the summed_up matrix
            all_solutions(j,k) = mat2cell(matrix,rows,columns); % Append weight*matrix to cell array
        end
        S = summed_up - A; % Calculate S
        objective_values(j,1) = sum(sum(S.^2)); % Calculate objective of corresponding S
    end
end


    
