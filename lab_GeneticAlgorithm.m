clear
clc
close all

%load data
tsp262= load("tsp262.mat").A;
tsp202= load("tsp202.mat").A;
gr17 = load("gr17.mat").A;
eil76 = load("eil76.mat").A;
ch150 = load("ch150.mat").A;
berlin52 = load("berlin52.mat").A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameters to change
graph = berlin52;
P = 100;
K = P/2;
l = 5;
crossover_prob = 0.8;
mutation_prob = 0.7;
maxit = 200;
[global_cost,mean_costs,min_costs] = genetic_algorithm(graph,P,K,l,crossover_prob,mutation_prob,maxit);
plot(1:size(mean_costs,2),mean_costs)
hold on
plot(1:size(min_costs,2),min_costs)
xlabel('iteration')
ylabel('cost')
legend('mean cost','min cost')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [global_cost,mean_costs,min_costs] = genetic_algorithm(graph,P,K,l,crossover_prob,mutation_prob,maxit)
   global_cost = Inf;
   mean_costs = [];
   min_costs = [];
   it = 0;
   tours = initial_tours(graph,P);
   all_costs = zeros(1,P);
    for i = 1:P
        all_costs(1,i) = tour_cost(tours(i,:),graph);
    end
    min_cost = min(all_costs);
    if min_cost < global_cost
        global_cost = min_cost;
    end
    mean_costs(end+1) = mean(all_costs);
    min_costs(end+1) = min_cost;
    while it < maxit
        [parents,parents_costs] = selection_mating(tours,all_costs,K);
        new_parents = [];
        while size(new_parents,1) < P
            [p1, p2] = selection_parents(parents,parents_costs,l);
            [o1, o2] = offspring_creation(p1,p2,crossover_prob);
            o1 = mutation(o1,mutation_prob);
            o2 = mutation(o2,mutation_prob);
            new_parents(end+1,:) = o1;
            new_parents(end+1,:) = o2;
        end
        tours = new_parents;
        for i = 1:P
            all_costs(1,i) = tour_cost(tours(i,:),graph);
        end
        min_cost = min(all_costs);
        if min_cost < global_cost
            global_cost = min_cost;
        end
        mean_costs(end+1) = mean(all_costs);
        min_costs(end+1) = min_cost;
        it = it+1;
        fprintf('iteration: %g, optimal cost: %g\n',it,global_cost)
    end
end



%{
indices = zeros(1,length(parents_costs));
all_p = 1:P;
for i = 1:length(parents_costs)
    element = parents_costs(i);
    index = all_costs == element;
    iteration = all_p(index);
    iteration = iteration(1);
    indices(1,i) = iteration;
end
scatter(indices,parents_costs)
xlabel('solution number')
ylabel('cost')
legend('all solutions','solutions chosen for mating pool')
%}


function tour = mutation(tour,mutation_prob)
    n = length(tour);
    u = rand;
    if u <= mutation_prob
        index = randi(n);
        if index == n
            first_element = tour(n);
            second_element = tour(1);
            tour(1) = first_element;
            tour(n) = second_element;
        else
            first_element = tour(index);
            second_element = tour(index + 1);
            tour(index) = second_element;
            tour(index + 1) = first_element;
        end
    end
end

function [offspring1, offspring2] = offspring_creation(parent1,parent2,crossover_prob)
    n = length(parent1);
    cutoff = round(n/2);
    all_vertices = 1:n;
    
    offspring1 = zeros(1,n);
    offspring1(1:cutoff) = parent1(1:cutoff);
    to_find = all_vertices(~ismember(all_vertices,parent1(1:cutoff)));
    offspring1(cutoff+1:end) = intersect(parent2,to_find,'stable');
    
    offspring2 = zeros(1,n);
    offspring2(1:cutoff) = parent2(1:cutoff);
    to_find2 = all_vertices(~ismember(all_vertices,parent2(1:cutoff)));
    offspring2(cutoff+1:end) = intersect(parent1,to_find2,'stable');
    
    u = rand;
    if u > crossover_prob
        offspring1 = parent1;
        offspring2 = parent2;
    end
end

function [parent1, parent2] = selection_parents(parents,parents_costs,l)
    possibilities = size(parents,1);
    mate_parents = randperm(possibilities,l);
    costs = parents_costs(mate_parents);
    [~, idx1] = max(costs);
    parent1 = parents(idx1,:);
    costs(idx1) = -Inf;
    [~, idx2] = max(costs);
    parent2 = parents(idx2,:);   
end

function [parents,parents_costs] = selection_mating(tours,all_costs,K)
    parents = [];
    parents_costs = [];
    P = length(all_costs);
    while size(parents,1) < K
        selected_parents = randi(P,1,K);
        all_parents = zeros(1,P);
        all_parents(selected_parents) = 1;
        selected_parents_indices = all_parents > 0;
        costs = all_costs(selected_parents_indices);
        mini = min(costs);
        [~,column] = find(all_costs == mini);
        column = column(1);
        parents(end+1,:) = tours(column,:);
        parents_costs(end+1) = mini;
    end


    %{
    n = size(tours,2);
    mini = min(all_costs);
    maxi = max(all_costs);
    prob = (all_costs - mini)./(maxi-mini);
    sel_prob = ones(1,length(all_costs)) - prob;
    parents = [];
    parents_costs = [];
  
    sol_used = [];
    while size(parents,1) < K
        parent = randi(size(tours,1));
        if sum(sol_used == parent) == 1
            continue
        else
            u = rand + 0.3;
            if u > sel_prob(parent)
                continue
            else
                sol_used(end+1) = parent;
                parents(end+1,:) = tours(parent,:);
                parents_costs(end+1) = all_costs(parent);
               
            end
        end     
    end
    %}
end

function tours = initial_tours(A,P)
    n = length(A);
    tours = zeros(P,n);
    for i = 1:P
        %{
        if i <= n
            tours(i,:) = initial_tour(A,i);
        else
        %}
            tours(i,:) = randperm(n);
        %end
    end
end

function cost = tour_cost(tour,A)
    %{
    %%%
    function calculates cost of a tour
    %%%
    input: tour of which you want to know the cost and distance matrix A
    %%%
    output: cost of the tour
    %}
    n = length(tour);
    cost = 0;
    for i=1:n-1
        cost = cost + A(tour(i),tour(i+1));
    end
    cost = cost + A(tour(n),tour(1));
end

function NN_tour = initial_tour(A,i)
    %{
    %%%
    function to create initial tour based on nearest neighbor heuristic
    %%%
    input: distance matrix A
    %%%
    output: nearest neighbor tour NN_tour
    %}
    n = size(A,1);
    NN_tour = [i];
    B = A;
    B(:,i) = 10000000;
    while length(NN_tour) < n
        minimum = min(B(NN_tour(end),:));
        new_vertex = find(B(NN_tour(end),:) == minimum);
        new_vertex = new_vertex(1);
        column_to_delete = new_vertex;
        NN_tour(end+1) = new_vertex;
        B(:,column_to_delete) = 10000000;
    end
end