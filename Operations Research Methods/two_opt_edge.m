clear
clc

tsp262= load("tsp262.mat").A;
tsp202= load("tsp202.mat").A;
gr17 = load("gr17.mat").A;
eil76 = load("eil76.mat").A;
ch150 = load("ch150.mat").A;
berlin52 = load("berlin52.mat").A;

graph = tsp262;

%% Exercise 1
tic
[tour_edge,cost_edge,it_edge] = iterative_TSP_edge(graph,-1); % -1 = trivial tour.
toc
disp("Iterative edge improvement gives opt. cost: "+num2str(cost_edge)+" for given TSP")

%% Exercise 2
tic
[tour_vertex,cost_vertex,it_vertex] = iterative_TSP_vertex(graph,-1); % -1 = trivial tour.
toc
disp("Iterative vertex improvement gives opt. cost: "+num2str(cost_vertex)+" for given TSP")

%% Exercise 3


%% Functions

function [tour,cost,it] = iterative_TSP_edge(A,starting_vertex)

    if starting_vertex == -1
        tour = 1:1:size(A,2); 
    else
        tour = NN_heuristic(A,starting_vertex); % nearest neighbor tour
    end

    cost = tour_cost(tour,A);
    stop = 0;
    it = 1;
    while stop == 0
        % find best tour j (=new_tour variable) in 2 edge exchange neighborhood of tour i (= tour variable)
        [new_tour,new_cost] = two_edge_exchange_neighborhood(A,tour,cost);
        if new_cost < cost
            cost = new_cost;
            tour = new_tour;
            it = it +1;
        else
            stop = 1;
        end
    end
end

function [tour,cost,it] = iterative_TSP_vertex(A,starting_vertex)
    if starting_vertex == -1
        tour = 1:1:size(A,2); 
    else
        tour = NN_heuristic(A,starting_vertex); % nearest neighbor tour
    end

    cost = tour_cost(tour,A);
    stop = 0;
    it = 1;
    while stop == 0
        % find best tour j (=new_tour variable) in 2 vertex exchange neighborhood of tour i (= tour variable)
        [new_tour,new_cost] = two_vertex_exchange_neighborhood(A,tour,cost);
        if new_cost < cost
            cost = new_cost;
            tour = new_tour;
            it = it +1;
        else
            stop = 1;
        end
    end
end

function [best_tour,best_cost] = two_edge_exchange_neighborhood(A,tour,cost)
    N = size(A,2);
    org_cost = cost; % original cost
    best_tour = tour; % output
    best_cost = cost; % output
    % we swap edge {a,b} with edge {c,d}
    iter = 0;
    for i=1:N-2
        a = tour(i);
        b = tour(i+1);
        for j=i+2:N
            c = tour(j);
            if j == N
                d = tour(1);
            else
                d = tour(j+1);
            end

            new_cost = cost - A(a,b) - A(c,d) + A(a,c) + A(b,d);
            new_tour = tour;
            if j==N
                reverse = flip(tour(1:i));
                new_tour(1:i) = reverse;
            else
                reverse = flip(tour(i+1:j));
                new_tour(i+1:j) = reverse;
            end

            if new_cost < best_cost
                best_cost = new_cost;
                best_tour = new_tour;
            end        
            iter = iter +1;
        end        
    end
end

function [best_tour,best_cost] = two_vertex_exchange_neighborhood(A,tour,cost)

    N = size(A,2);
    best_tour = tour; % store the currently best improvement in the 2 vertex neighborhood of tour
    best_cost = cost; % store the currently best cost corresponding to best_tour

    for i=1:N-1
        % determine predecessor and succesor of vertex i
        if i == 1
            i_minus = tour(N);
            i_plus = tour(i+1);
        else
            i_plus = tour(i+1);
            i_minus = tour(i-1);
        end

        for j=i+1:N
            % determine predecessor and succesor of vertex j
            if j == N
                j_minus = tour(j-1);
                j_plus = tour(1);
            else
                j_minus = tour(j-1);
                j_plus = tour(j+1);
            end

            % vertex i swaps with vertex j
            new_tour = tour;
            new_tour(i) = tour(j);
            new_tour(j) = tour(i);
            
            % new cost = old cost - outgoing edges + ingoing edges
            if i_plus == tour(j) % i is directly behind j in tour
                new_cost = cost - A(i_minus,tour(i)) - A(tour(j),j_plus) + A(i_minus,tour(j)) + A(tour(i),j_plus);
            elseif j_plus == tour(i) % j is directly behind i in tour
                new_cost = cost - A(j_minus,tour(j)) - A(tour(i),i_plus) + A(j_minus,tour(i)) + A(tour(j),i_plus);
            else % i and j are not directly to each other
                new_cost = cost - A(i_minus,tour(i)) - A(tour(i),i_plus) -A(j_minus,tour(j)) - A(tour(j),j_plus) + A(i_minus,tour(j))  +A(tour(j),i_plus) + A(j_minus,tour(i))+ A(tour(i),j_plus);
            end

            if new_cost < best_cost
                best_tour = new_tour;
                best_cost = new_cost;
            end
        end
    end
end

function tour = NN_heuristic(A, starting_vertex)
    N = size(A,2);
    visited(1:N) = 0;
    visited(starting_vertex) = 1;
    tour(1:N) = 0; %preallocation
    tour(1) = starting_vertex;
    for i=1:N-1
        row = A(tour(i), :);
        row(~~visited) = Inf;
        [~, index] = min(row);
        visited(index) = 1;
        tour(i+1) = index;
    end
end

function cost = tour_cost(tour,A)
    N = size(tour,2);
    cost = 0;
    for i=1:N-1
        cost = cost + A(tour(i),tour(i+1));
    end
    cost = cost + A(N,1);
end