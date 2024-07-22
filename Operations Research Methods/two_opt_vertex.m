%function two_opt_vertex(A)
tic
    clear
    clc
    A = load('gr17.mat');
    A = A.A;
    n = size(A,1);
    all_vertices = 1:n;
    initial_vertex = 1;
    NN_tour = [1];
    B = A;
    B(:,1) = 10000000;
    while length(NN_tour) < n
        minimum = min(B(NN_tour(end),:));
        new_vertex = find(B(NN_tour(end),:) == minimum);
        new_vertex = new_vertex(1);
        column_to_delete = new_vertex;
        NN_tour(end+1) = new_vertex;
        B(:,column_to_delete) = 10000000;
    end
    
    for j = 1:n-1
        begin_point = NN_tour(j);
        end_point = NN_tour(j+1);
        x(begin_point,end_point) = 1;
        x(end_point,begin_point) = 1;
        if j == n-1
            begin_point = NN_tour(n);
            end_point = NN_tour(1);
            x(begin_point,end_point) = 1;
            x(end_point,begin_point) = 1;
        end
    end
    
    current_cost = sum(sum(x.*A))/2;
    fprintf('current cost is %g\n',current_cost)
    possibilities = nchoosek(n,2);
    alpha = 1;
    while alpha == 1
        all_tours = zeros(possibilities,n);
        k = 1;
        for i = 1:n-1
            for j = i+1:n
                new_tour = NN_tour;
                first_element = NN_tour(i);
                second_element = NN_tour(j);
                new_tour(i) = second_element;
                new_tour(j) = first_element;
                all_tours(k,:) = new_tour;
                k = k+1;
            end
        end

        costs = [];
        matrices = {};

        for m = 1:possibilities
            row = all_tours(m,:);
            x = zeros(n);
            for j = 1:n-1
                begin_point = row(j);
                end_point = row(j+1);
                x(begin_point,end_point) = 1;
                x(end_point,begin_point) = 1;
                if j == n-1
                    begin_point = row(n);
                    end_point = row(1);
                    x(begin_point,end_point) = 1;
                    x(end_point,begin_point) = 1;
                end
            end
            cost = sum(sum(x.*A))/2;
            costs(end+1) = cost;
            matrices{end+1} = x;
        end
        min_cost = min(costs);
        if min_cost < current_cost
            current_cost = min_cost;
            index = find(costs == min_cost);
            index = index(1);
            new_matrix = cell2mat(matrices(index));
            x = new_matrix;

            new_tour_opt = [1];
            new_matrix_copy = new_matrix;
            new_matrix_copy(:,1) = 0;
            while length(new_tour_opt) < n
                positions = find(new_matrix_copy(new_tour_opt(end),:) == 1);
                position = positions(1);
                new_tour_opt(end+1) = position;
                new_matrix_copy(:,position) = 0;
            end
            NN_tour = new_tour_opt;
            fprintf('a new better tour has been found: %g\n',min_cost)
        else
           alpha = 0; 
           fprintf('no improvement was possible, current tour is optimal tour with cost %g\n',current_cost) 
        end
    end
toc
%end