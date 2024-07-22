clear
clc
close all

GCP1 = load('GCP1.mat').A;
GCP2 = load('GCP2.mat').A;
GCP3 = load('GCP3.mat').A;
GCP4 = load('GCP4.mat').A;
GCP5 = load('GCP5.mat').A;
GCP6 = load('GCP6.mat').A;
GCP7 = load('GCP7.mat').A;

graph = GCP6;

tic
[coloring,chromatic_numbers] = GCP_tabu(graph);
fprintf('chromatic number is %g\n',chromatic_numbers(end))
toc


function [coloring,all_chromatic_numbers] = GCP_tabu(A)
    %A = graph;
    %A = [0 1 1 0;1 0 0 1;1 0 0 1;0 1 1 0];
    n = size(A,2);
    all_vertices = 1:n;
    ini_coloring = initial_coloring(A);
    max_color = max(ini_coloring);
    %max color indicates the maximum number of colors needed to find a
    %feasible solution and is therefore an upper bound for the problem
    %which can then be used to construct the first (infeasible) coloring
    coloring = ones(1,n);%round(1+(max_color-1)*rand(1,n));
    coloring = coloring(1:n);
    all_chromatic_numbers = [max_color];
    iter = 0;
    gamma_matrix = zeros(n,max_color);
    tabu_matrix = zeros(n,max_color);
    for i = 1:n
        adjacent_vertices = A(i,:) > 0;
        colors = coloring(adjacent_vertices);
        for j = 1:length(colors)
            element = colors(1,j);
            gamma_matrix(i,element) = gamma_matrix(i,element) + 1;
        end
    end
    
   [number_conflicting_edges,tabu_increaser,conflicting_vertices] = conflicting_edges(gamma_matrix,coloring,A,iter);
   all_numbers_conflicting_edges = [number_conflicting_edges];
    while number_conflicting_edges ~= 0 && iter < 100000
        tabu_check = 0;
        iter = iter + 1;
        if mod(iter,10000) == 0
            fprintf('iteration: %g\n',iter)
        end
        delta_matrix = gamma_matrix;
        non_conflicting_vertices = all_vertices(~ismember(all_vertices,conflicting_vertices));
        delta_matrix(non_conflicting_vertices,:) = Inf;
        for i = 1:length(conflicting_vertices)
            current_vertex = conflicting_vertices(i);
            current_color = coloring(current_vertex);
            substraction_term = gamma_matrix(current_vertex,current_color);
            delta_matrix(current_vertex,:) = delta_matrix(current_vertex,:) - substraction_term;
        end
        
        alpha = 1;
        while alpha == 1
            min_delta = min(min(delta_matrix));
            [row,column] = find(delta_matrix == min_delta);
            row = row(1);
            column = column(1);
            if tabu_matrix(row,column) > iter
                delta_matrix(row,column) = Inf;
                new_coloring = coloring;
                new_coloring(row) = column;
                row_indices = A(row,:) > 0;
                all_vertices = 1:n;
                rows_to_change = all_vertices(row_indices);
                tabu_color = coloring(row);
                gamma_matrix_check = gamma_matrix;
                gamma_matrix_check(rows_to_change,tabu_color) = gamma_matrix_check(rows_to_change,tabu_color) - 1;
                gamma_matrix_check(rows_to_change,column) = gamma_matrix_check(rows_to_change,column) + 1;
                [number_conflicting_edges_check,tabu_increaser_check,conflicting_vertices_check] = conflicting_edges(gamma_matrix_check,coloring_copy,A,iter);
                if number_conflicting_edges_check == 0
                    tabu_check = 1;
                    alpha = 0;
                    gamma_matrix_new = gamma_matrix_check;
                    number_conflicting_edges_new = number_conflicting_edges_check;
                end
            else
                alpha = 0;
            end
        end
        if tabu_check == 0
            row_indices = A(row,:) > 0;
            all_vertices = 1:n;
            rows_to_change = all_vertices(row_indices);
            tabu_color = coloring(row);
            tabu_matrix(row,tabu_color) = tabu_matrix(row,tabu_color) + tabu_increaser;
            gamma_matrix_new = gamma_matrix;
            gamma_matrix_new(rows_to_change,tabu_color) = gamma_matrix_new(rows_to_change,tabu_color) - 1;
            gamma_matrix_new(rows_to_change,column) = gamma_matrix_new(rows_to_change,column) + 1;
            new_coloring = coloring;
            new_coloring(row) = column;
            [number_conflicting_edges_new,~,conflicting_vertices_new] = conflicting_edges(gamma_matrix_new,new_coloring,A,iter);
        end
        
        all_numbers_conflicting_edges(end+1) = number_conflicting_edges_new;
        if number_conflicting_edges_new < number_conflicting_edges
            number_conflicting_edges = number_conflicting_edges_new;
            gamma_matrix = gamma_matrix_new;
            conflicting_vertices = conflicting_vertices_new;
            coloring = new_coloring;
            all_chromatic_numbers(end+1) = max(coloring);
        else
            all_chromatic_numbers(end+1) = all_chromatic_numbers(end);
        end
    end
 
 yyaxis left
 grid on
 all_chromatic_numbers = all_chromatic_numbers(1:iter+1);
 plot(1:iter+1,all_chromatic_numbers)
 hold on
 plot(1:iter+1,repmat(max_color,1,iter+1),'--r')
 xlabel('iterations')
 ylabel('chromatic number')
 hold on
 yyaxis right
 plot(1:iter+1,all_numbers_conflicting_edges)
 ylabel('number of conflicting edges')
 title('chromatic number and number of conflicting edges')
 legend('chromatic number','upper bound for max number of colors','number of conflicting edges')
 

 end
function coloring = initial_coloring(A)
    n = size(A,2);
    coloring = zeros(1,n);
    coloring(1) = 1;
    for i = 2:n
        adjacent_vertices = A(i,:) > 0;
        colors = coloring(adjacent_vertices);
        new_color = max(colors)+1;
        coloring(1,i) = new_color;
    end
end

function [number_conflicting_edges,tabu_increaser,conflicting_vertices] = conflicting_edges(gamma_matrix,coloring,A,iter)
    n = size(gamma_matrix,1);
    max_color = size(gamma_matrix,2);
    F_c = gamma_matrix(sub2ind([n, max_color],1:n,coloring));
    conflicting_vertices_indices = F_c > 0;
    tabu_increaser = iter + 5 + 0.1*sum(conflicting_vertices_indices);
    all_vertices = 1:n;
    conflicting_vertices = all_vertices(conflicting_vertices_indices);
    conflicting_colors = coloring(1,conflicting_vertices);
    conflicted_vertices_colors = [conflicting_vertices;conflicting_colors];
    max_conflicting_color = max(conflicting_colors);
    number_conflicting_edges = 0;
    for i = 1:max_conflicting_color
        [rows, columns] = find(conflicted_vertices_colors(2,:) == i);
        if isempty(columns) ~= 1
            vertices_to_inspect = conflicted_vertices_colors(1,columns);
            number_conflicting_edges = number_conflicting_edges + sum(sum(A(vertices_to_inspect,vertices_to_inspect)));
        else
            continue
        end
    end
    
end
