clear
clc
close all

rng(189)

% Initialization:
M = load('inst3.mat').M;
K_max = 12; % Max # elements on right hand side
W_max = 150; % Max sum of weights
T = 10; % Stopping criterion: time

P = 1500;
l = 800;
mating = 500;
l2 = 50;

if exist('P','var') == 1
    [best_cost,global_opt,all_k,all_w,B,S,k,w] = ARCDP2(M,K_max,W_max,T,P,l,mating,l2);
else
    [best_cost,global_opt,all_k,all_w,B,S,k,w] = ARCDP2(M,K_max,W_max,T);
end

function [best_cost,global_opt,all_k,all_w,B,S,k,w] = ARCDP2(A,K_max,W_max,T,P,l,mating,l2)
    sameit = 20; % stopping criterion: stop if after sameit iterations global opt has not changed
    mutation_prob = 0.2;
    maxit = 5000;

    if nargin == 4
        P = max(200,2000-1/5*(size(A,1)*size(A,2))); % # initial solutions
        l = round(0.5*P); % size tournament mating pool
        mating = round(P/3); % size mating pool
        l2 = round(0.1 * mating); % size tournament parent selection
    end

    tstart = tic;
    best_cost = zeros(1,maxit);
    all_k = zeros(1,maxit);
    all_w = zeros(1,maxit);
    global_opt = Inf;
    it = 1;
    totaltime = 0;
    iteration_same = 0;
    stop = 1;
    [initial_solution,objective_values] = initial_solutions(A,K_max,W_max,P);
    minimum = min(objective_values);
    best_cost(1,it) = minimum;
    if minimum < global_opt
        global_opt = minimum;
        idx = find(objective_values == global_opt);
        idx = idx(1);
        B = initial_solution(idx,:);
        all_k(1,it) = size(B,2);
        all_w(1,it) = sum(cell2mat(cellfun(@(x) max(max(x)),B,'UniformOutput',false)));
    end
    tend = toc(tstart);
    totaltime = totaltime + tend;
    fprintf('----------- %g -----------\n',it)
    fprintf('best cost : %g\n',minimum)
    fprintf('global opt: %g\n',global_opt)
    fprintf('total time: %g\n',totaltime)
    fprintf('--------------------------\n')
    fprintf('\n')
    while it < maxit && global_opt ~= 0 && totaltime < T && stop == 1 %print == 1 is stopping criterion for sameit: if it > sameit => print = 0 so the new iteration will not be printed in the command window
        tstart = tic;
        it = it+1;
        [mating_pool,mating_pool_obj_values] = mating_pool_creator(objective_values,mating,l);
        parents = parents_selection(mating_pool,mating_pool_obj_values,P,l2);
        [children,objective_values] = offspring_creation(A,initial_solution,P,parents,K_max,W_max,mutation_prob);
        initial_solution = children;
        minimum = min(objective_values);
        best_cost(1,it) = minimum;
        if minimum < global_opt
            iteration_same = 0;
            global_opt = minimum;
            idx = find(objective_values == global_opt);
            idx = idx(1);
            B = initial_solution(idx,:);
            B = B(1,~cellfun(@isempty,B));
            all_k(1,it) = size(B,2);
            all_w(1,it) = sum(cell2mat(cellfun(@(x) max(max(x)),B,'UniformOutput',false)));
        else
            iteration_same = iteration_same + 1;
            idx = find(objective_values == minimum);
            idx = idx(1);
            C = initial_solution(idx,:);
            C = C(1,~cellfun(@isempty,C));
            all_k(1,it) = size(C,2);
            all_w(1,it) = sum(cell2mat(cellfun(@(x) max(max(x)),C,'UniformOutput',false)));
            if iteration_same > sameit
                stop = 0;
            end
        end
        tend = toc(tstart);
        totaltime = totaltime + tend;
        fprintf('----------- %g -----------\n',it)
        fprintf('best cost : %g\n',minimum)
        fprintf('global opt: %g\n',global_opt)
        fprintf('total time: %g\n',totaltime)
        fprintf('--------------------------\n')
        fprintf('\n')
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
    subplot(2,1,1)
    plot(1:it,best_cost(1:it))
    subplot(2,1,2)
    plot(1:it, all_w(1:it))
end

function [children,objective_values] = offspring_creation(A,initial_solution,P,parents,K_max,W_max,mutation_prob)
    num_matrices_vector = sum(~cellfun(@isempty,initial_solution),2); % Find # of matrices of each solution
    rows = size(initial_solution{1,1},1);
    columns = size(initial_solution{1,1},2);
    children = cell(P,K_max);
    objective_values = zeros(P,1);
    iter = 1;
    for row = 1:round(P/2) 
        sum_child1 = zeros(rows,columns);
        sum_child2 = zeros(rows,columns);
        
        couple = parents(row,:);
        couple = sort(couple);
        parent1 = couple(1);
        parent2 = couple(2);
        num_matrices = num_matrices_vector(couple,:); % Find # matrices of both parents

        if num_matrices(1) > num_matrices(2)
            [parent2,parent1] = deal(parent1,parent2);
        end
        
        % Create cell of the parents, the matrices of the parents and the
        % weights of the parents
        parent_cell = [initial_solution(parent1,:);initial_solution(parent2,:)];
        parent_matrices_cell = cellfun(@(x) double(x>0),parent_cell,'UniformOutput', false);
        parent_weights = cellfun(@(x) max(max(x)),parent_cell,'UniformOutput',false);
        
        num_matrices = num_matrices_vector([parent1,parent2],:);
        max_matrices = max(num_matrices);
        min_matrices = min(num_matrices);

        
        weight_vector1 = cell2mat(parent_weights(1,:));
        weight_vector2 = cell2mat(parent_weights(2,:));
        
        weight_vector1(end+1:max_matrices) = 0; 

        % Initialize cell for offspring:
        offspring = cell(2,K_max);

        % Creating offspring:
        j = 1;
        for i = 1:max_matrices
            matrix1 = parent_matrices_cell{1,j};
            matrix2 = parent_matrices_cell{2,i};
            offspring_matrix1 = matrix1;
            offspring_matrix2 = matrix2;

            % Determine the weight of both matrices
            childweight1 = weight_vector1(i);
            childweight2 = weight_vector2(i);
            
            
            % Swap a random row of the matrices to create two new offspring matrices
            row_tbs = randperm(rows,round(rows/2));
            offspring_matrix1(row_tbs,:) = matrix2(row_tbs,:);
            offspring_matrix2(row_tbs,:) = matrix1(row_tbs,:);
            
            % Apply random mutation: swap a random row within a child
            u = rand;
            if u <= mutation_prob
                row_swap = randperm(rows,2);
                % Swap rows of first child
                stored_rows = offspring_matrix1(row_swap,:);
                offspring_matrix1(row_swap(1),:) = stored_rows(2,:);
                offspring_matrix1(row_swap(2),:) = stored_rows(1,:);
                % Swap rows of second child:
                stored_rows = offspring_matrix2(row_swap,:);
                offspring_matrix2(row_swap(1),:) = stored_rows(2,:);
                offspring_matrix2(row_swap(2),:) = stored_rows(1,:);
            end

            % Add found matrices to the child 3D matrices
            offspring{1,i} = childweight1.*offspring_matrix1;
            offspring{2,i} = childweight2.*offspring_matrix2;
           
            % Add found matrices to a sum 
            sum_child1 = sum_child1 + childweight1.*offspring_matrix1;
            sum_child2 = sum_child2 + childweight2.*offspring_matrix2;
            
            if j < min_matrices % Makes sure that if we reach the end of parent 1's
                                % matrices we go back to the first matrix
                j = j+1;
            else
                j = 1;
            end
        end
        children(2*iter-1,:) = offspring(1,:);
        children(2*iter,:) = offspring(2,:);
        
        % Calculate objective values of found offspring
        objective_values(2*iter-1,1) = sum(sum((sum_child1 - A).^2));
        objective_values(2*iter,1) = sum(sum((sum_child2 - A).^2));
        iter = iter + 1;

    end
end

function parents = parents_selection(mating_pool,mating_pool_obj_values,P,l2)
    % Selects parent couples from the mating pool and stores these pairs in
    % a matrix. Takes as input variables a row vector as mating pool, # of 
    % initial solutions and size of the tournament.
    mating = length(mating_pool);
    parents = zeros(round(P/2),2);
    
    i = 1;
    while i <= round(P/2)
        indices = randperm(mating,l2); % Choose l2 distinct random indices between 1 and mating
        values = mating_pool_obj_values(indices); %Find corresponding values of selected indices
        [~,row] = min(values);
        mating_pool_obj_values(indices(row)) = Inf; %Select index with the minimal value and set its value to +Inf
        indices2 = randperm(mating,l2);
        values2 = mating_pool_obj_values(indices2);
        [~,row2] = min(values2);
        % Append the indices of both parents to the matrix "parents"
        parents(i,1) = mating_pool(indices(row));
        parents(i,2) = mating_pool(indices2(row2));
        i = i+1;
    end
 end



 function [mating_pool,mating_pool_obj_values] = mating_pool_creator(objective_values,mating,l)
    % Creates a mating pool of size M by tournament selection. Takes as
    % input variables the objective values of the solutions,the size of the
    % mating pool and the size of the tournament.
    P = length(objective_values);
    mating_pool = zeros(1,mating);
    objective_values_copy = objective_values;
    
    i = 1;
    mating_pool_obj_values = zeros(1,mating);
    while i <= mating
        indices = randperm(P,l); % Select l district random indices between 1 and P
        values = objective_values_copy(indices); %Find corresponding objective values of selected indices
        [mini,row] = min(values);
        mating_pool_obj_values(1,i) = mini;
        mating_pool(1,i) = indices(row); % Add the index of the solution with the minimal value in the mating pool 
        objective_values_copy(indices(row)) = Inf; % Set value of selected solution to Inf to make sure it doesn't get added in the mating pool again
        i = i+1;
    end
 end


function [all_solutions,objective_values] = initial_solutions(A,K_max,W_max,P)
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
            max_weight = weight_left - weight - matrices_left + 1 + 1;
            weight_left = weight_left - weight;
            matrices_left = matrices_left - 1;

            matrix = matrix.*weight; % Multiply matrix by its weight to easily store the weight
            summed_up = summed_up + matrix; % Add new found weight*matrix to the summed_up matrix
            all_solutions(j,k) = mat2cell(matrix,rows,columns); % Append weight*matrix to cell array
        end
        S = summed_up - A; % Calculate S
        objective_values(j,1) = sum(sum(S.^2)); % Calculate objective of corresponding S
    end
end
