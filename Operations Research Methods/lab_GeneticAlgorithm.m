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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameters to change
graph = gr17;
P = 50;
K = 25;
l = 5;
crossover_prob = 0.8;
mutation_prob = 0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
tours = initial_tours(graph,P);
all_costs = zeros(1,P);
for i = 1:P
    all_costs(1,i) = tour_cost(tours(i,:),graph);
end
%scatter(1:P,all_costs)
%hold on
[parents,parents_costs] = selection_mating(tours,all_costs,K);
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
[p1, p2] = selection_parents(parents,parents_costs,l);
toc

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
    n = size(tours,2);
    mini = min(all_costs);
    maxi = max(all_costs);
    prob = (all_costs - mini)./(maxi-mini);
    sel_prob = ones(1,length(all_costs)) - prob;
    parents = zeros(K,n);
    parents_costs = zeros(1,K);
    counter = 1;
    sol_used = [];
    while counter <= K
        parent = randi(length(tours));
        if sum(sol_used == parent) == 1
            continue
        else
            u = rand;
            if u > sel_prob(parent)
                continue
            else
                sol_used(counter) = parent;
                parents(counter,:) = tours(parent,:);
                parents_costs(1,counter) = all_costs(parent);
                counter = counter + 1;
            end
        end     
    end     
end

function tours = initial_tours(A,P)
    n = length(A);
    tours = zeros(P,n);
    for i = 1:P
        tours(i,:) = randperm(n);
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