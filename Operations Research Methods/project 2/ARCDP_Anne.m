clear
clc
close all

% Initialization:
A = [1 2 3; 4 5 6]; % Matrix
K_max = 10; % Max # elements on right hand side
W_max = 10; % Max sum of weights
P = 50; % # initial solutions
l = 5; % size tournament mating pool
M = 20; % size mating pool
l2 = 3; % size tournament parent selection


[initial_solution,objective_values] = initial_solutions(A,K_max,W_max,P)
mating_pool = mating_pool_creater(objective_values,M,l)
parents = parents_selection(mating_pool,P,l2)

%%% CROSSOVER SECTION %%% (change into function later)
num_matrices_vector = sum(~cellfun(@isempty,initial_solution),2); % Find # of matrices of each solution
rows = size(initial_solution{1,1},1);
columns = size(initial_solution{1,1},2);

% When we change this into a function we can give the row in parents as
% input variable
couple = parents(1,:);
couple = sort(couple);
parent1 = couple(1);
parent2 = couple(2);
num_matrices = num_matrices_vector(couple,:); % Find # matrices of both parents
    
matrix3D_parent1 = zeros(rows,columns,num_matrices(1));
weights_parent1 = zeros(1,num_matrices(1));
for j = 1:num_matrices(1)
    matrix = initial_solution{parent1,j};
    matrix3D_parent1(:,:,j) = double(matrix>0); % Add parent 1's matrix to 3D matrix
    weights_parent1(1,j) = max(matrix(matrix>=0)); % Add parent 1's weight to weight vector
end

matrix3D_parent2 = zeros(rows,columns,num_matrices(2));
weights_parent2 = zeros(1,num_matrices(1));
for j = 1:num_matrices(2)
    matrix = initial_solution{parent2,j};
    matrix3D_parent2(:,:,j) = double(matrix>0); % Add parent 1's matrix to 3D matrix
    weights_parent2(1,j) = max(matrix(matrix>=0)); % Add parent 1's weight to weight vector
end

max_matrices = max(num_matrices);
min_matrices = min(num_matrices);

% Initialize 3D matrices for offspring:
child1 = zeros(rows,columns,max_matrices);
child2 = zeros(rows,columns,max_matrices);

% Creating offspring:
j = 1;
for i = 1:max_matrices
    matrix1 = matrix3D_parent1(:,:,j);
    matrix2 = matrix3D_parent2(:,:,i);
    offspring_matrix1 = matrix1;
    offspring_matrix2 = matrix2;
    
    % Swap a random row of the matrices to create to new offspring matrices
    row_tbs = randi(rows);
    offspring_matrix1(row_tbs,:) = matrix2(row_tbs,:);
    offspring_matrix2(row_tbs,:) = matrix1(row_tbs,:);
    
    % Add found matrices to the child 3D matrices
    child1(:,:,i) = offspring_matrix1;
    child2(:,:,i) = offspring_matrix2;
        
    if j < min_matrices % Makes sure that if we reach the end of parent 1's
                        % matrices we go back to the first matrix
        j = j+1;
    else
        j = 1;
    end
end


 function parents = parents_selection(mating_pool,P,l2)
    % Selects parent couples from the mating pool and stores these pairs in
    % a matrix. Takes as input variables a row vector as mating pool, # of 
    % initial solutions and size of the tournament.
    M = length(mating_pool);
    parents = zeros(P/2,2);
    
    i = 1;
    while i <= P/2
        copy_mating_pool = mating_pool;
        indices = randperm(M,l2); % Choose l2 distinct random indices between 1 and M
        values = mating_pool(indices); %Find corresponding values of selected indices
        [~,row] = min(values);
        copy_mating_pool(indices(row)) = Inf; %Select index with the minimal value and set its value to +Inf
        indices2 = randperm(M,l2);
        values2 = copy_mating_pool(indices2);
        [~,row2] = min(values2);
        % Append the indices of both parents to the matrix "parents"
        parents(i,1) = indices(row);
        parents(i,2) = indices2(row2);
        i = i+1;
    end
 end



 function mating_pool = mating_pool_creater(objective_values,M,l)
    % Creates a mating pool of size M by tournament selection. Takes as
    % input variables the objective values of the solutions,the size of the
    % mating pool and the size of the tournament.
    P = length(objective_values);
    mating_pool = zeros(1,M);
    
    i = 1;
    while i <= M
        indices = randperm(P,l); % Select l district random indices between 1 and P
        values = objective_values(indices); %Find corresponding objective values of selected indices
        [~,row] = min(values);
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

    
