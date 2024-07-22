function [optimal_value,used_budget,new_budget,upgrade_matrix_shsr,upgrade_matrix_hsr,cost_upgrade_matrix_shsr,cost_upgrade_matrix_hsr] = one_shortest_path_shortened(budget,alpha,nodes,rhs_till_now,obj,sparse_constraint_matrix,number_routes,number_edges,cost_shsr,cost_hsr)
    %%%output:
    %optimal_value: maximum_shortest_path
    %used_budget: budget used
    %new_budget: optimized budget for given budget and max shortest path
    %upgrade_matrix_shsr: edges upgraded to shsr in initial round
    %upgrade_matrix_hsr: edges upgraded to hsr in initial round
    %cost_upgrade_matrix_shsr: edges upgraded to shsr in cost minimization
    %cost_upgrade_matrix_hsr: edges upgraded to hsr in cost minimization

    
    addpath("C:\gurobi1001\win64\matlab")

    rhs_fifth_constraint = budget;
    rhs = [rhs_till_now; rhs_fifth_constraint];
    model = [];
    model.obj = obj;
    model.A = sparse_constraint_matrix;
    model.rhs = rhs;
    model.sense = [repmat('>',number_routes,1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<'];
    model.modelsense = 'min';
    model.vtype = ['C';repmat('B',2*number_edges,1)];
    params.outputflag = 0;
    
    solution = gurobi(model,params);
    optimal_value = solution.objval;
    
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
    
    new_obj = [0 cost_shsr cost_hsr];
    new_constraint = [1 zeros(1,2*number_edges)];
    new_rhs = optimal_value*(1+alpha);

    cost_model = [];
    cost_model.obj = new_obj;
    cost_model.A = [sparse_constraint_matrix; new_constraint];
    cost_model.rhs = [rhs; new_rhs];
    cost_model.sense = [repmat('>',number_routes,1); repmat('<',number_edges,1); repmat('<',number_edges,1); repmat('<', number_edges,1); '<';'<'];
    cost_model.modelsense = 'min';
    cost_model.vtype = ['C';repmat('B',2*number_edges,1)];
    params.outputflag = 0;

    cost_solution = gurobi(cost_model,params);
    if strcmp(cost_solution.status,'INFEASIBLE') == 0
        new_budget = cost_solution.objval;
        new_length_maximum_shortest_path = cost_solution.x(1);

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
        new_budget = 0;
        cost_upgrade_matrix_shsr = zeros(49);
        cost_upgrade_matrix_hsr = zeros(49);
    end
end

