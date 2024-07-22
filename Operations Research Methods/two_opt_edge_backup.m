%function [new_tour, NN_tour] = two_opt_edge(A)
    %%%computes local 2-opt minimum by swapping edges
    %input: A - cost_matrix
    %output: new_tour - optimal tour found
    %        NN_tour - initial_tour
    %        function will also print values to keep track of what is
    %        happening
    %--> this function also finds initial tour, namely nearest neighbour
    %tour starting in vertex 1
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
    
    possibilities = n*(n-3);
    start_vertex = NN_tour(1);
    last_vertex = NN_tour(end);
    no_first_vertices = [NN_tour(1) NN_tour(2) last_vertex];
    possible_first_vertices = setdiff(all_vertices,no_first_vertices);
    [~,indices,~] = intersect(NN_tour,possible_first_vertices);
    second_vertices = NN_tour(indices + 1);
    second_edges = [possible_first_vertices; second_vertices]';
    first_edge = repmat([NN_tour(1) NN_tour(2)],length(possible_first_vertices),1);
    second_edge = second_edges;
    
    no_first_vertices = [NN_tour(end) NN_tour(1) NN_tour(end-1)];
    possible_first_vertices = setdiff(all_vertices,no_first_vertices);
    [~,indices,~] = intersect(NN_tour,possible_first_vertices);
    second_vertices = NN_tour(indices + 1);
    second_edges = [possible_first_vertices;second_vertices]';
    first_edges = repmat([NN_tour(end) NN_tour(1)],length(possible_first_vertices),1);
    
    first_edge = [first_edge;first_edges];
    second_edge = [second_edge;second_edges];
    
    NN_tour_helper = [NN_tour 1];
    k = 2;
    while size(first_edge,1) < possibilities
        first_vertex_loop = NN_tour(k);
        no_first_vertices = [NN_tour(k) NN_tour(k+1) NN_tour(k-1)];
        possible_first_vertices = setdiff(all_vertices,no_first_vertices);
        [~,indices,~] = intersect(NN_tour,possible_first_vertices);
        second_vertices = NN_tour_helper(indices + 1);
        second_edges = [possible_first_vertices;second_vertices]';
        first_edges = repmat([NN_tour(k) NN_tour(k+1)],length(possible_first_vertices),1);

        first_edge = [first_edge;first_edges];
        second_edge = [second_edge;second_edges];
        k = k+1;

    end
    
    x = zeros(n);
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
    fprintf('current cost is: %g\n', current_cost)
    
    costs = [];
    matrices = {};
    
    for i = 1:size(first_edge,1)
        y = x;
        first_deleted_edge = first_edge(i,:);
        second_deleted_edge = second_edge(i,:);
        vertex1_1 = first_deleted_edge(1);
        vertex1_2 = first_deleted_edge(2);
        index1_1 = find(NN_tour == vertex1_1);
        
        vertex2_1 = second_deleted_edge(1);
        vertex2_2 = second_deleted_edge(2);
        index2_1 = find(NN_tour == vertex2_1);
        
        y(vertex1_1,vertex1_2) = 0;
        y(vertex1_2,vertex1_1) = 0;
        
        y(vertex2_1,vertex2_2) = 0;
        y(vertex2_2,vertex2_1) = 0;
        
        y(vertex1_1,vertex2_1) = 1;
        y(vertex2_1,vertex1_1) = 1;
        
        y(vertex1_2,vertex2_2) = 1;
        y(vertex2_2,vertex1_2) = 1;
        
        cost = sum(sum(y.*A))/2;
        costs(end+1) = cost;
        matrices{end+1} = y;
    end
    
    min_cost = min(costs);
    if min_cost < current_cost
        index = find(costs == min_cost);
        index = index(1);
        new_matrix = cell2mat(matrices(index));
        x = new_matrix;

        new_tour = [1];
        new_matrix_copy = new_matrix;
        new_matrix_copy(:,1) = 0;
        while length(new_tour) < n
            positions = find(new_matrix_copy(new_tour(end),:) == 1);
            position = positions(1);
            new_tour(end+1) = position;
            new_matrix_copy(:,position) = 0;
        end
        fprintf('a new better tour has been found: %g\n',min_cost)
        new_tour;
    else
       fprintf('no improvement was possible, initial tour is optimal tour\n') 
    end
 
    %%%vanaf hier tweede functie
    while min_cost < current_cost
        start_vertex = new_tour(1);
        last_vertex = new_tour(end);
        no_first_vertices = [new_tour(1) new_tour(2) last_vertex];
        possible_first_vertices = setdiff(all_vertices,no_first_vertices);
        [~,indices,~] = intersect(new_tour,possible_first_vertices);
        second_vertices = new_tour(indices + 1);
        second_edges = [possible_first_vertices; second_vertices]';
        first_edge = repmat([new_tour(1) new_tour(2)],length(possible_first_vertices),1);
        second_edge = second_edges;

        no_first_vertices = [new_tour(end) new_tour(1) new_tour(end-1)];
        possible_first_vertices = setdiff(all_vertices,no_first_vertices);
        [~,indices,~] = intersect(new_tour,possible_first_vertices);
        second_vertices = new_tour(indices + 1);
        second_edges = [possible_first_vertices;second_vertices]';
        first_edges = repmat([new_tour(end) new_tour(1)],length(possible_first_vertices),1);

        first_edge = [first_edge;first_edges];
        second_edge = [second_edge;second_edges];

        new_tour_helper = [new_tour 1];
        k = 2;
        while size(first_edge,1) < possibilities
            first_vertex_loop = new_tour(k);
            no_first_vertices = [new_tour(k) new_tour(k+1) new_tour(k-1)];
            possible_first_vertices = setdiff(all_vertices,no_first_vertices);
            [~,indices,~] = intersect(new_tour,possible_first_vertices);
            second_vertices = new_tour_helper(indices + 1);
            second_edges = [possible_first_vertices;second_vertices]';
            first_edges = repmat([new_tour(k) new_tour(k+1)],length(possible_first_vertices),1);

            first_edge = [first_edge;first_edges];
            second_edge = [second_edge;second_edges];
            k = k+1;

        end
        
        current_cost = min_cost;

        costs = [];
        matrices = {};

        for i = 1:size(first_edge,1)
            y = x;
            first_deleted_edge = first_edge(i,:);
            second_deleted_edge = second_edge(i,:);
            vertex1_1 = first_deleted_edge(1);
            vertex1_2 = first_deleted_edge(2);
            index1_1 = find(NN_tour == vertex1_1);

            vertex2_1 = second_deleted_edge(1);
            vertex2_2 = second_deleted_edge(2);
            index2_1 = find(NN_tour == vertex2_1);

            y(vertex1_1,vertex1_2) = 0;
            y(vertex1_2,vertex1_1) = 0;

            y(vertex2_1,vertex2_2) = 0;
            y(vertex2_2,vertex2_1) = 0;

            y(vertex1_1,vertex2_1) = 1;
            y(vertex2_1,vertex1_1) = 1;

            y(vertex1_2,vertex2_2) = 1;
            y(vertex2_2,vertex1_2) = 1;

            cost = sum(sum(y.*A))/2;
            costs(end+1) = cost;
            matrices{end+1} = y;
        end

        min_cost = min(costs);
        if min_cost < current_cost
            index = find(costs == min_cost);
            index = index(1);
            new_matrix = cell2mat(matrices(index));
            x = new_matrix;

            new_tour = [1];
            new_matrix_copy = new_matrix;
            new_matrix_copy(:,1) = 0;
            while length(new_tour) < n
                positions = find(new_matrix_copy(new_tour(end),:) == 1);
                position = positions(1);
                new_tour(end+1) = position;
                new_matrix_copy(:,position) = 0;
            end
            fprintf('a new better tour has been found: %g\n', min_cost)
            new_tour;
        else
           fprintf('no improvement was possible, current tour is optimal tour with value %g\n',current_cost) 
        end

    end 
%end