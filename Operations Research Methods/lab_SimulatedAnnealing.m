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
graph = ch150;
T = 110;
alpha = 0.999;
I = 310;
maxit = 1000000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
%best_cost = Inf;
%best_T = 0;
%best_I = 0;
%for T=100:10:500
    %for I = 100:10:500
        [tour,cost,all_costs,all_temperatures,all_acceptance_prob,global_minimum] = Simulated_Annealing(graph,T,alpha,I,maxit);
        %fprintf('T = %g, I = %g, optimum = %g\n',T,I,cost)
        %{if cost < best_cost
            %best_cost = cost;
            %best_T = T;
            %best_I = I;
        %end
    %end
%end
 
fprintf('global minimum is %g\n',global_minimum)
%toc

%create plots
subplot(2,2,1)
x = 1:1:length(all_costs);
plot(x,all_costs)
xlabel('iteration')
ylabel('cost')
title('cost per iteration')

subplot(2,2,2)
x2 = 1:1:length(all_temperatures);
plot(x2,all_temperatures);
xlabel('iteration')
ylabel('temperature')
title('temperature per iteration')

subplot(2,2,3)
x3 = 1:1:length(all_acceptance_prob);
scatter(x3,all_acceptance_prob)
xlabel('iteration')
ylabel('acceptance probability')
title('acceptance probability per iteration')

function NN_tour = initial_tour(A)
    %{
    %%%
    function to create initial tour based on nearest neighbor heuristic
    %%%
    input: distance matrix A
    %%%
    output: nearest neighbor tour NN_tour
    %}
    n = size(A,1);
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

function new_tour = neighborhood_searcher(tour)
    %{
    %%%
    function searches for new tour by swapping two edges
    %%%
    input: tour for which you want to swap two edges
    %%%
    output: new tour in which two random edges are swapped
    %}
    new_tour = tour;
    n = length(tour);
    
    
    i = round(1 + (n-1) * rand(1,1));
    if i == 1
        vertices_left = [3:n-1];
    elseif i == 2
        vertices_left = [4:n];
    elseif i == n-1
        vertices_left = [1:n-3];
    elseif i == n
        vertices_left = [2:n-2];
    else
        vertices_left = [1:i-2,i+2:n];
    end
    m = length(vertices_left);
    index = randperm(m,1);
    j = vertices_left(index);
    i_star = min(i,j);
    j_star = max(i,j);
    if j_star == n
        reverse = flip(tour(1:i_star));
        new_tour(1:i_star) = reverse;
    else
        reverse = flip(tour(i_star+1:j_star));
        new_tour(i_star+1:j_star) = reverse;
    end
    
end


function [tour,cost,all_costs,all_temperatures,all_acceptance_prob,global_minimum] = Simulated_Annealing(A,T,alpha,I,maxit)
    %{
    %%%
    function calculates optimal tour given user inputted parameter
    %%%
    input: A: distance matrix
           T: initial temperature level
           alpha: T = alpha*T after I iterations, usually 0.9 or 0.99
           I: number of iterations after which temperature is decreased
           maxit: maximum number of iterations algorithm can run
    %%%
    output: tour: most recent optimal tour
            cost: most recent optimal cost
            all_costs,all_temperatures&all_accptance_prob: arrays used in
            plotting after function is finished
            global_minimum: best found solution, if eg. 500 has been found
            earlier and last optimum found is 600, global_minimum will
            return 500 whereas cost will return 600
    %}
    i = 0;
    global_minimum = Inf;
    tour = initial_tour(A);
    old_cost = tour_cost(tour,A);
    all_costs = [old_cost];
    all_temperatures = [T];
    all_acceptance_prob = [];
    it = 0;
    while it < maxit
        new_tour = neighborhood_searcher(tour);
        new_cost = tour_cost(new_tour,A);
        if new_cost < old_cost
            acceptance_prob = 1;
            if new_cost < global_minimum
                global_minimum = new_cost;
            end
        else
            acceptance_prob = min(1,exp((old_cost - new_cost)/T));
        end
        all_acceptance_prob(end+1) = acceptance_prob;
        u = rand;
        if u <= acceptance_prob
            tour = new_tour;
            old_cost = tour_cost(tour,A);
        end
        all_costs(end+1) = old_cost;
        if i == I
            T = alpha*T;
            i = 0;
        else
            i = i+1;
        end
        all_temperatures(end+1) = T;
        it = it+1;       
    end
    cost = tour_cost(tour,A);
end