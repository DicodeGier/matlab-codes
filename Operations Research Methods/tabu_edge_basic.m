%this version is probably incorrect due to the neighborhood chosen but
%arrives at the (almost) correct solutions

clear
clc
close all

tsp262= load("tsp262.mat").A;
tsp202= load("tsp202.mat").A;
gr17 = load("gr17.mat").A;
eil76 = load("eil76.mat").A;
ch150 = load("ch150.mat").A;
berlin52 = load("berlin52.mat").A;

graph = berlin52;

tic
[tour_edge,cost_edge,it_edge,tabu_matrix] = iterative_TSP_edge(graph,1); % -1 = trivial tour.
toc

%for ch150: tabu_increaser = 50 and step_increaser = 4
%otherwise step_increaser = 3 and tabu_increaser = 30

function [tour,cost,it,tabu_matrix] = iterative_TSP_edge(A,starting_vertex)
    n = size(A,2);    
    tabu_increaser = 30;
    times_increasing = 0;
    step_increaser = 3;
    if starting_vertex == -1
        tour = 1:1:size(A,2); 
    else
        tour = NN_heuristic(A,starting_vertex); % nearest neighbor tour
    end

    cost = tour_cost(tour,A);
    global_optimum = 2*cost;
    fprintf('current cost is: %g\n',cost)
    tabu_matrix = zeros(n);
    stop = 0;
    it = 1;
    all_costs = [cost];
    while it < 2000
        % find best tour j (=new_tour variable) in 2 edge exchange neighborhood of tour i (= tour variable)
        [new_tour,new_cost,best_a,best_b,best_c,best_d,times_increasing] = two_edge_exchange_neighborhood(A,tour,cost,tabu_matrix,global_optimum,times_increasing,step_increaser);
        all_costs(end+1) = new_cost;
        %if new_cost < cost
            if new_cost < global_optimum
                global_optimum = new_cost;
                elements = ones(n);
                to_delete = triu(elements,1);
                elements_to_delete = to_delete > 0;
                tabu_matrix(elements_to_delete) = 0;
            
            %lower triangular part
            elseif best_a ~= 0
                tabu_matrix(max([best_a,best_b]),min([best_a,best_b])) = tabu_matrix(max([best_a,best_b]),min([best_a,best_b])) + 1;
                tabu_matrix(max([best_c,best_d]),min([best_c,best_d])) = tabu_matrix(max([best_c,best_d]),min([best_c,best_d])) + 1;
                    
                %upper triangular part
                tabu_matrix(min([best_a,best_b]),max([best_a,best_b])) = tabu_matrix(min([best_a,best_b]),max([best_a,best_b])) + tabu_increaser;
                tabu_matrix(min([best_c,best_d]),max([best_c,best_d])) = tabu_matrix(min([best_c,best_d]),max([best_c,best_d])) + tabu_increaser;
                %tabu_matrix(min([best_a,best_c]),max([best_a,best_c])) = tabu_matrix(min([best_a,best_c]),max([best_a,best_c])) + tabu_increaser;
                %tabu_matrix(min([best_b,best_d]),max([best_b,best_d])) = tabu_matrix(min([best_b,best_d]),max([best_b,best_d])) + tabu_increaser;

                if it > 1
                    elements_larger_zero = tabu_matrix > 0;
                    elements_to_reduce = triu(elements_larger_zero,1);
                    tabu_matrix(elements_to_reduce) = tabu_matrix(elements_to_reduce) - 1;
                end
          
                if new_cost < cost
                    fprintf('a better tour has been found: %g, iteration: %g\n', new_cost,it)
                    %check = tour_cost(new_tour,A)
                else
                    fprintf('cost is increasing to get over the hill: %g,iteration: %g\n',new_cost,it)
                    %check = tour_cost(new_tour,A)
                end
                cost = new_cost;
                tour = new_tour;
                it = it +1;
            else
                fprintf('no more possibilities to try, global optimum is %g\n',global_optimum)
                all_costs = all_costs(1:it);
                plot(1:it,all_costs)
                cost = new_cost;
                tour = new_tour;
                return
            end
        %{
        else
            stop = 1;
            
   
        end
            %}
    end
     fprintf('maximum number of iterations run, global optimum is %g\n',global_optimum)
     all_costs = all_costs(1:it);
     plot(1:it,all_costs)
     xlabel('iterations')
     ylabel('costs')
end


function [best_tour,best_cost,best_a,best_b,best_c,best_d,times_increasing] = two_edge_exchange_neighborhood(A,tour,cost,tabu_matrix,global_optimum,times_increasing,step_increaser)
    N = size(A,2);
    best_cost = cost; % original cost
    best_tour = tour; % output    
    if times_increasing > 0
      times_increasing = times_increasing - 1;
      [best_tour,best_cost,best_a,best_b,best_c,best_d] = two_edge_exchange_neigborhood_increase(A,best_tour,best_cost,tabu_matrix);  
    else
    best_a = 0;
    best_b = 0;
    best_c = 0;
    best_d = 0;
    % we swap edge {a,b} with edge {c,d}
    iter = 0;
    for i=1:N-2
        a = tour(i);
        b = tour(i+1);
        for j=i+2:N
            tabu_edge = 0;
            c = tour(j);
            if j == N
                d = tour(1);
            else
                d = tour(j+1);
            end
            if tabu_matrix(min([a,b]),max([a,b])) > 0 | tabu_matrix(min([c,d]),max([c,d])) > 0
                tabu_edge = 1;
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
            if tabu_edge == 1
                if new_cost < best_cost - 20 && new_cost ~= cost
                    best_a = a;
                    best_b = b;
                    best_c = c;
                    best_d = d;
                    best_cost = new_cost;
                    best_tour = new_tour;
                end      
            elseif new_cost < best_cost && new_cost ~= cost
                best_a = a;
                best_b = b;
                best_c = c;
                best_d = d;
                best_cost = new_cost;
                best_tour = new_tour;
            end
        end
        iter = iter +1;
    end
    if best_a == 0
        [best_tour,best_cost,best_a,best_b,best_c,best_d] = two_edge_exchange_neigborhood_increase(A,best_tour,best_cost,tabu_matrix);
        times_increasing = times_increasing + step_increaser;
    end
    end
 end

function [best_tour,best_cost_upper,best_a,best_b,best_c,best_d] = two_edge_exchange_neigborhood_increase(A,tour,cost,tabu_matrix)
    N = size(A,2);
    best_cost = cost; % original cost
    best_cost_upper = Inf;
    best_tour = tour; % output
    best_a = 0;
    best_b = 0;
    best_c = 0;
    best_d = 0;
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
            if tabu_matrix(min([a,b]),max([a,b])) > 0 | tabu_matrix(min([c,d]),max([c,d])) > 0
                continue
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
            if new_cost > best_cost && new_cost < best_cost_upper
                best_a = a;
                best_b = b;
                best_c = c;
                best_d = d;
                best_cost_upper = new_cost;
                best_tour = new_tour;
            end
        end
        iter = iter +1;
    end 
end

function tour = NN_heuristic(A, starting_vertex)
    tour = [starting_vertex];
    A(:,starting_vertex) = 10000000;
    n = size(A,2);
    while length(tour) < n
       minimum = min(A(tour(end),:));
       new_vertex = find(A(tour(end),:) == minimum);
       new_vertex = new_vertex(1);
       column_to_delete = new_vertex;
       tour(end+1) = new_vertex;
       A(:,column_to_delete) = 10000000;
    end
end

function cost = tour_cost(tour,A)
    N = size(tour,2);
    cost = 0;
    for i=1:N-1
        cost = cost + A(tour(i),tour(i+1));
    end
    cost = cost + A(tour(N),tour(1));
end
