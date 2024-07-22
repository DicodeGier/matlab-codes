function [cost_difference, cost_upgrade_matrix_shsr, cost_upgrade_matrix_hsr] = one_shortest_path(budget, alpha)
    addpath("C:\gurobi1001\win64\matlab")

    %%%load data and construct parameters
    data = importfile('data network complete.xlsx', 'tijden (uren)');
    hours = data.data;
    clear data

    data = importfile('data network complete.xlsx', 'afstanden (kilometers)');
    distances = data.data;
    clear data

    %check if data is complete
    for i = 1:size(hours,1)
        for j = i+1:size(hours,1)
            if hours(i,j) > 0
                if distances(i,j) > 0
                    continue
                else
                    fprintf('error at %g, %g\n',i,j)
                end
            end
        end
    end

    indices = triu(ones(size(hours,1)),1);
    indices = indices > 0;

    data = importfile('data network complete.xlsx', 'u_shsr');
    u_shsr = data.data;
    u_shsr(isnan(u_shsr)) = 0;
    clear data
    %u_shsr_matrix = u_shsr;
    u_shsr(~indices) = Inf;
    u_shsr = u_shsr';
    u_shsr = u_shsr(:);
    u_shsr(u_shsr == Inf) = [];
    u_shsr = u_shsr';

    data = importfile('data network complete.xlsx', 'u_hsr');
    u_hsr = data.data;
    u_hsr(isnan(u_hsr)) = 0;
    clear data
    %u_hsr_matrix = u_hsr;
    u_hsr(~indices) = Inf;
    u_hsr = u_hsr';
    u_hsr = u_hsr(:);
    u_hsr(u_hsr == Inf) = [];
    u_hsr = u_hsr';

    %distances_matrix = distances;
    distances(~indices) = Inf;
    distances = distances';
    distances = distances(:);
    distances(distances == Inf) = [];
    distances = distances';

    t = hours;
    t(~indices) = Inf;
    t = t';
    t = t(:);
    t(t == Inf) = [];
    t = t';

    shsr = u_shsr .* distances./250;
    hsr = u_hsr .* distances./300;

    cost_shsr = 16.*distances.*u_shsr;
    cost_hsr = 20.*distances.*u_hsr;

    u_shsr_0 = u_shsr == 0;
    u_shsr_1 = u_shsr == 1;
    u_shsr_copy = u_shsr;
    u_shsr_copy(u_shsr_0) = 1;
    u_shsr_copy(u_shsr_1) = 0;
    shsr = shsr + u_shsr_copy .* t;

    u_hsr_0 = u_hsr == 0;
    u_hsr_1 = u_hsr == 1;
    u_hsr_copy = u_hsr;
    u_hsr_copy(u_hsr_0) = 1;
    u_hsr_copy(u_hsr_1) = 0;
    hsr = hsr + u_hsr_copy .* t;

    shsr_matrix = zeros(size(hours,1));
    shsr_copy = shsr;
    for i = 1:size(hours,1)
       for j = 1:size(hours,1)
           if j <= i
               shsr_matrix(i,j) = 0;
           else
               shsr_matrix(i,j) = shsr_copy(1);
               shsr_copy(1) = [];
           end
       end
    end

    hsr_matrix = zeros(size(hours,1));
    hsr_copy = hsr;
    for i = 1:size(hours,1)
       for j = 1:size(hours,1)
           if j <= i
               hsr_matrix(i,j) = 0;
           else
               hsr_matrix(i,j) = hsr_copy(1);
               hsr_copy(1) = [];
           end
       end
    end

    red_shsr = t - shsr;
    red_hsr = t - hsr;

    red_shsr_matrix = zeros(size(hours,1));
    red_shsr_copy = red_shsr;
    for i = 1:size(hours,1)
       for j = 1:size(hours,1)
           if j <= i
               red_shsr_matrix(i,j) = 0;
           else
               red_shsr_matrix(i,j) = red_shsr_copy(1);
               red_shsr_copy(1) = [];
           end
       end
    end

    red_hsr_matrix = zeros(size(hours,1));
    red_hsr_copy = red_hsr;
    for i = 1:size(hours,1)
       for j = 1:size(hours,1)
           if j <= i
               red_hsr_matrix(i,j) = 0;
           else
               red_hsr_matrix(i,j) = red_hsr_copy(1);
               red_hsr_copy(1) = [];
           end
       end
    end


    u_hsr = u_hsr';
    u_shsr = u_shsr';

    %%%construct the shortest paths
    nodes = 49;
    real_nodes = 35; %nodes between which the shortest path is calculated (i.e. not the artificial nodes)
    possibilities = real_nodes*(real_nodes - 1)/2; %maximum number of shortest paths
    result = cell(possibilities,6);

    counter = 1;
    for i = 1:nodes
        %exclude artificial nodes
        if i == 5 || i == 6 || i == 7 || i == 9 || i == 16 || i == 17 || i == 18 || i == 23 || i == 24 || i == 29 || i == 30 || i == 31 || i == 32 || i == 44
            continue
        end
        for j = i+1:nodes
            if j == 5 || j == 6 || j == 7 || j == 9 || j == 16 || j == 17 || j == 18 || j == 23 || j == 24 || j == 29 || j == 30 || j == 31 || j == 32 || j == 44
                continue
            end
            %calculate shortest path
            [cost, route] = dijkstra(hours,i,j);

            %store required information in a cell
            route = flip(route,2);
            result{counter,1} = i; %starting node
            result{counter,2} = j; %ending node
            result{counter,3} = route; %traversed route
            result{counter,4} = cost; %cost of route
            matrix = zeros(nodes);
            for k = 1:length(route)-1
               element_1 = route(k);
               element_2 = route(k+1);
               row_node = min(element_1,element_2);
               column_node = max(element_1,element_2);
               matrix(row_node,column_node) = 1;
            end
            result{counter,5} = matrix; %matrix in which edges are in the upper triangular part
            matrix(~indices) = Inf;
            matrix = matrix';
            matrix = matrix(:);
            matrix(matrix==Inf) = [];
            matrix = matrix';
            result{counter,6} = matrix; %vector of edges stretched along the second dimension, i.e. x = [x11 x12 x13 ...] so "the rows first"
            counter = counter + 1;

        end
    end

    %hist([result{:,4}])
    max_length = max([result{:,4}]);

    %%%GEHARDCODED: DIT NOG DYNAMISCHER MAKEN!!!!%%%
    %exclude routes that are less than the mean travel time
    number_routes = size(result,1);
    number_edges = length(result{1,6});
    x_st = zeros(number_routes,number_edges);
    for k = 1:number_routes
        x_st(k,:) = result{k,6};
    end

    %objective function
    z_obj = 1;
    y_1_obj = zeros(1,number_edges);
    y_2_obj = zeros(1,number_edges);

    obj = [z_obj y_1_obj y_2_obj];

    %first constraint: shortest path longer than any path after upgrading
    y_1_first_constraint = red_shsr.*x_st;
    y_2_first_constraint = red_hsr.*x_st;
    z_first_constraint = ones(number_routes,1);
    first_constraint_matrix = [z_first_constraint y_1_first_constraint y_2_first_constraint];
    rhs_first_constraint = t*transpose(x_st);

    %second constraint: edge can only be upgraded once
    y_1_second_constraint = eye(number_edges);
    y_2_second_constraint = eye(number_edges);
    z_second_constraint = zeros(number_edges,1);
    second_constraint_matrix = [z_second_constraint y_1_second_constraint y_2_second_constraint];
    rhs_second_constraint = ones(number_edges,1);

    %third constraint: only upgradeable edges are upgraded to shsr
    z_third_constraint = zeros(number_edges,1);
    y_1_third_constraint = eye(number_edges);
    y_2_third_constraint = zeros(number_edges);
    third_constraint_matrix = [z_third_constraint y_1_third_constraint y_2_third_constraint];
    rhs_third_constraint = u_shsr;

    %fourth constraint: only upgradeable edges are upgraded to hsr
    z_fourth_constraint = zeros(number_edges,1);
    y_1_fourth_constraint = zeros(number_edges);
    y_2_fourth_constraint = eye(number_edges);
    fourth_constraint_matrix = [z_fourth_constraint y_1_fourth_constraint y_2_fourth_constraint];
    rhs_fourth_constraint = u_hsr;


    %fifth constraint: cost within budget
    z_fifth_constraint = 0;
    y_1_fifth_constraint = cost_shsr;
    y_2_fifth_constraint = cost_hsr;
    fifth_constraint_matrix = [z_fifth_constraint y_1_fifth_constraint y_2_fifth_constraint];
    rhs_fifth_constraint = budget;

    constraint_matrix = [first_constraint_matrix; second_constraint_matrix; third_constraint_matrix; fourth_constraint_matrix; fifth_constraint_matrix];
    rhs = [transpose(rhs_first_constraint); rhs_second_constraint; rhs_third_constraint; rhs_fourth_constraint; rhs_fifth_constraint];

    %gurobi model
    model = [];
    model.obj = obj;
    model.A = sparse(constraint_matrix);
    model.rhs = rhs;
    model.sense = [repmat('>',number_routes,1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<'];
    model.modelsense = 'min';
    model.vtype = ['C';repmat('B',2*number_edges,1)];
    params.outputflag = 0;

    solution = gurobi(model,params);
    optimal_value = solution.objval;

    %create matrix with edges that are upgraded for better visualization
    upgrades = solution.x;
    upgrades(1) = [];
    upgrades_shsr = upgrades(1:number_edges);
    upgrades_hsr = upgrades(number_edges+1:end);
    
    used_budget = cost_shsr*upgrades_shsr + cost_hsr*upgrades_hsr;

    upgrades_shsr_copy = upgrades_shsr;
    upgrades_hsr_copy = upgrades_hsr;
    upgrade_matrix_hsr = zeros(nodes);
    row_counter = 0;
    for i = 1:nodes
       row_counter = row_counter + 1;
       for j = 1:nodes
          if j + row_counter > nodes
              break
          end
          upgrade_matrix_hsr(i,j+row_counter) = upgrades_hsr_copy(1);
          upgrades_hsr_copy(1) = [];
       end
    end

    upgrade_matrix_shsr = zeros(nodes);
    row_counter = 0;
    for i = 1:nodes
       row_counter = row_counter + 1;
       for j = 1:nodes
          if j + row_counter > nodes
              break
          end
          upgrade_matrix_shsr(i,j+row_counter) = upgrades_shsr_copy(1);
          upgrades_shsr_copy(1) = [];
       end
    end
    
    %now minimize the cost given the optimal value
    new_obj = [0 cost_shsr cost_hsr];
    new_constraint = [1 zeros(1,2*number_edges)];
    new_rhs = optimal_value*(1+alpha);

    cost_model = [];
    cost_model.obj = new_obj;
    cost_model.A = sparse([constraint_matrix; new_constraint]);
    cost_model.rhs = [rhs; new_rhs];
    cost_model.sense = [repmat('>',number_routes,1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<';'<'];
    cost_model.modelsense = 'min';
    cost_model.vtype = ['C';repmat('B',2*number_edges,1)];
    params.outputflag = 0;

    cost_solution = gurobi(cost_model,params);
    if strcmp(cost_solution.status,'INFEASIBLE') == 0
        optimal_cost = cost_solution.objval;
        new_length_maximum_shortest_path = cost_solution.x(1);
        percentual_increase = (new_length_maximum_shortest_path - optimal_value)/optimal_value * 100;

        cost_upgrades = cost_solution.x;
        cost_upgrades(1) = [];
        cost_upgrades_shsr = cost_upgrades(1:number_edges);
        cost_upgrades_hsr = cost_upgrades(number_edges+1:end);

        cost_upgrades_shsr_copy = cost_upgrades_shsr;
        cost_upgrades_hsr_copy = cost_upgrades_hsr;
        cost_upgrade_matrix_hsr = zeros(nodes);
        row_counter = 0;
        for i = 1:nodes
           row_counter = row_counter + 1;
           for j = 1:nodes
              if j + row_counter > nodes
                  break
              end
              cost_upgrade_matrix_hsr(i,j+row_counter) = cost_upgrades_hsr_copy(1);
              cost_upgrades_hsr_copy(1) = [];
           end
        end

        cost_upgrade_matrix_shsr = zeros(nodes);
        row_counter = 0;
        for i = 1:nodes
           row_counter = row_counter + 1;
           for j = 1:nodes
              if j + row_counter > nodes
                  break
              end
              cost_upgrade_matrix_shsr(i,j+row_counter) = cost_upgrades_shsr_copy(1);
              cost_upgrades_shsr_copy(1) = [];
           end
        end
    else
        fprintf('error with budget of %g and alpha of %g\n',budget,alpha)
        optimal_cost = 0;
        cost_upgrade_matrix_shsr = zeros(49);
        cost_upgrade_matrix_hsr = zeros(49);
    end
    cost_difference = used_budget - optimal_cost;

end

