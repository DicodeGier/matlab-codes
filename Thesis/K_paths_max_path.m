function [optimal_value, used_budget, upgrade_matrix_shsr, upgrade_matrix_hsr, all_routes_selected] = K_paths_max_path(budget,rhs_till_now,obj,sparse_constraint_matrix,number_routes,number_edges,cost_shsr,cost_hsr,all_K,number_binary_variables)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IMPORTANT:
    %Notice that K is not an input in the function, this function only
    %works when the constraints, given K, have been prespecified
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rhs_fifth_constraint = budget;
    rhs = [rhs_till_now; rhs_fifth_constraint];
    model = [];
    model.obj = obj;
    model.A = sparse_constraint_matrix;
    model.rhs = rhs;
    model.sense = [repmat('>',number_routes,1); repmat('=',length(all_K),1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<'];
    model.modelsense = 'min';
    model.vtype = ['C';repmat('B',2*number_edges,1);repmat('B',number_binary_variables,1)];
    params.outputflag = 0;

    solution = gurobi(model,params);
    optimal_value = solution.objval;

    upgrades = solution.x;
    upgrades(1) = [];
    upgrades_shsr = upgrades(1:number_edges);
    upgrades_hsr = upgrades(number_edges+1:2*number_edges);

    used_budget = cost_shsr*upgrades_shsr + cost_hsr*upgrades_hsr;

    upgrades_shsr_copy = upgrades_shsr;
    upgrades_hsr_copy = upgrades_hsr;
    nodes = 49;
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
    all_routes_selected = zeros(length(all_K),1);
    
    counter = 1;
    routes_selected = solution.x(2*number_edges+2:end);
    routes_selected = round(routes_selected,2);
    for z = 1:length(all_K)
        routes = all_K(z);
        routes_selected_subset = routes_selected(1:routes);
        [row,~] = find(routes_selected_subset == 0);
        routes_selected(1:routes) = [];
        all_routes_selected(counter,1) = row;
        counter = counter + 1;
    end
    
end
