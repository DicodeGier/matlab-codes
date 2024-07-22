clear
clc
close all

rng(1)

% Initialization:
M = load('inst3.mat').M; %[3,6,4;2,1,5];%load('inst3.mat').M;
K_max = 16; % Max # elements on right hand side
W_max = 200; % Max sum of weights
P = 1000; % # initial solutions
l = 50; % size tournament mating pool
M_size = 500; % size mating pool
l2 = 25; % size tournament parent selection
maxit = 50; % stopping criterion: maximum of 1000 iterations
sameit = 3000; %s stopping criterion: stop if after sameit iterations global opt has not changed
mutation_prob = 0.2;

tic
[best_cost,global_opt] = ARCDP2(M,K_max,W_max,P,M_size,l,l2,maxit,sameit,mutation_prob);
toc

fprintf('global_optimum is %g\n',global_opt)
plot(1:maxit,best_cost)

function [best_cost,global_opt] = ARCDP2(A,K_max,W_max,P,M_size,l,l2,maxit,sameit,mutation_prob)
    best_cost = zeros(1,maxit);
    global_opt = Inf;
    it = 1;
    iteration_same = 0;
    [initial_solution,objective_values] = initial_solutions(A,K_max,W_max,P);
    minimum = min(objective_values);
    best_cost(1,it) = minimum;
    if minimum < global_opt
        global_opt = minimum;
    end
    fprintf('iteration: %g\n',it)
    while it < maxit && global_opt ~= 0
        it = it+1;
        [mating_pool,mating_pool_obj_values] = mating_pool_creater(objective_values,M_size,l);
        parents = parents_selection(mating_pool,mating_pool_obj_values,P,l2);
        [children,objective_values] = offspring_creation(A,initial_solution,P,parents,K_max,W_max,mutation_prob);
        initial_solution = children;
        minimum = min(objective_values);
        best_cost(1,it) = minimum;
        if minimum < global_opt
            global_opt = minimum;
        else
            iteration_same = iteration_same + 1;
            if iteration_same > sameit
                return
            end
        end
        fprintf('iteration: %g\n',it)
    end
end

function [children,objective_values] = offspring_creation(A,initial_solution,P,parents,K_max,W_max,mutation_prob)
    num_matrices_vector = sum(~cellfun(@isempty,initial_solution),2); % Find # of matrices of each solution
    rows = size(initial_solution{1,1},1);
    columns = size(initial_solution{1,1},2);
    children_3d = cell(P,1);
    objective_values = zeros(P,1);
    iter = 1;
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
        
        matrix3D_parent1 = zeros(rows,columns,num_matrices(1));
        weights_parent1 = zeros(1,num_matrices(1));
        for j = 1:num_matrices(1)
            matrix = initial_solution{parent1,j};
            matrix3D_parent1(:,:,j) = double(matrix>0); % Add parent 1's matrix to 3D matrix
            weights_parent1(1,j) = max(matrix(matrix>=0)); % Add parent 1's weight to weight vector
        end

        matrix3D_parent2 = zeros(rows,columns,num_matrices(2));
        weights_parent2 = zeros(1,num_matrices(2));
        for j = 1:num_matrices(2)
            matrix = initial_solution{parent2,j};
            matrix3D_parent2(:,:,j) = double(matrix>0); % Add parent 1's matrix to 3D matrix
            weights_parent2(1,j) = max(matrix(matrix>=0)); % Add parent 1's weight to weight vector
        end

        max_matrices = max(num_matrices);
        min_matrices = min(num_matrices);
        
        matrices_to_use = round((max_matrices + min_matrices)/2);

        % Initialize 3D matrices for offspring:
        child1 = zeros(rows,columns,matrices_to_use);
        child2 = zeros(rows,columns,matrices_to_use);

        % Creating offspring:
        j = 1;
        for i = 1:matrices_to_use
            matrix1 = matrix3D_parent1(:,:,j);
            matrix2 = matrix3D_parent2(:,:,i);
            offspring_matrix1 = matrix1;
            offspring_matrix2 = matrix2;

            % determine the weight of both matrices
            weight1 = weights_parent1(j);
            weight2 = weights_parent2(i);
            childweight = (weight1 + weight2)/2;
            if mod(childweight,1) ~= 0
                childweight = round(childweight) - 1;
            end
            
            u = rand;
            if u <= mutation_prob
                if childweight == W_max - max_matrices
                    childweight = childweight - 1;
                else
                    childweight = childweight + 1;
                end
            end

            % Swap a random row of the matrices to create to new offspring matrices
            all_rows = 1:rows;
            n = round(rows/2);
            residual1 = sum((A - matrix1).^2,2);
            rows1 = zeros(1,n);
            for i = 1:n
                [~,idx] = min(residual1);
                rows1(i) = idx;
                residual1(idx) = Inf;
            end          
            not_rows1 = all_rows(~ismember(all_rows,rows1));
            
            residual2 = sort(sum((A - matrix2).^2,2));
            rows2 = zeros(1,n);
            for i = 1:n
                [~,idx] = min(residual2);
                rows2(i) = idx;
                residual2(idx) = Inf;
            end 
            not_rows2 = all_rows(~ismember(all_rows,rows2));
            
            offspring_matrix1(not_rows1,:) = matrix2(rows2,:);
            offspring_matrix2(not_rows2,:) = matrix1(rows1,:);

            % Add found matrices to the child 3D matrices
            child1(:,:,i) = childweight.*offspring_matrix1;
            child2(:,:,i) = childweight.*offspring_matrix2;

            if j < min_matrices % Makes sure that if we reach the end of parent 1's
                                % matrices we go back to the first matrix
                j = j+1;
            else
                j = 1;
            end
        end
        children_3d{2*iter-1} = child1;
        objective_values(2*iter-1,1) = sum(sum((sum(child1,3) - A).^2));
        children_3d{2*iter} = child2;
        objective_values(2*iter,1) = sum(sum((sum(child2,3) - A).^2));
        iter = iter + 1;

    end
    children = cell(P,K_max);
    for i = 1:P
        dimensions = size(cell2mat(children_3d(i,1)),3);
        element = cell2mat(children_3d(i,1));
        for j = 1:dimensions
            children{i,j} = element(:,:,j);
        end
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

    
