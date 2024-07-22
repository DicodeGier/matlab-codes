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

graph = GCP7;
[coloring,chromatic_number] = tabu_GCP(graph);
fprintf('chromatic number is %g\n',chromatic_number)

function [coloring,chromatic_number] = tabu_GCP(A)
    [coloring,k_star] = initial_coloring(A);
    check = 1;
    k = k_star;
    chromatic_number = k;
    fprintf('k is %g\n',k)
    while check == 1
        k = k-1;
        fprintf('trying k = %g\n',k)
        [coloring,check] = tabu_GCP_check(A,k);
        if check == 1
            chromatic_number = max(coloring);
        end
    end
end

function [coloring,check] = tabu_GCP_check(A,k)
    n = size(A,1);
    all_vertices = 1:n;
    coloring = round(1+(k-1)*rand(1,n));
    iter = 0;
    gamma_matrix = zeros(n,k);
    tabu_matrix = zeros(n,k);
    for i = 1:n
        adjacent_vertices = A(i,:) > 0;
        colors = coloring(adjacent_vertices);
        for j = 1:length(colors)
            element = colors(1,j);
            gamma_matrix(i,element) = gamma_matrix(i,element) + 1;
        end
    end
    
    [number_conflicting_edges,tabu_increaser,conflicting_vertices] = conflicting_edges(gamma_matrix,coloring,A,iter);
    while number_conflicting_edges ~= 0 && iter < 50000
        tabu_check = 0;
        iter = iter + 1;
        %{
        if mod(iter,10000) == 0
            fprintf('iteration: %g\n',iter)
        end
        %}
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
                [number_conflicting_edges_check,tabu_increaser_check,conflicting_vertices_check] = conflicting_edges(gamma_matrix_check,new_coloring,A,iter);
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
        
        if number_conflicting_edges_new < number_conflicting_edges
            number_conflicting_edges = number_conflicting_edges_new;
            gamma_matrix = gamma_matrix_new;
            conflicting_vertices = conflicting_vertices_new;
            coloring = new_coloring;
        end
    end
    if number_conflicting_edges == 0
        check = 1;
    else
        check = 0;
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

function [coloring,k_star] = initial_coloring(A)
    n = size(A,2);
    coloring = zeros(1,n);
    coloring(1) = 1;
    for i = 2:n
        adjacent_vertices = A(i,:) > 0;
        colors = coloring(adjacent_vertices);
        new_color = max(colors)+1;
        coloring(1,i) = new_color;
    end
    k_star = max(coloring);
end