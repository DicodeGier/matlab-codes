clear
clc
addpath("C:\gurobi1001\win64\matlab")

% Set parameters
truck_price = 300;
truck_size = 100;
penalty_coefficient = 0.2;
days = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'};

% Dict with demand info
demand = containers.Map({'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday'}, ...
                        [25, 75, 125, 225, 450, 250, 150]);

% Construct an empty optimization model
MIQP = gurobi.Model();

% Add integer variables to the model
delivery = MIQP.addVars(days, 'vtype', 'I', 'lb', 0, 'name', 'delivery');
trucks = MIQP.addVars(days, 'vtype', 'I', 'lb', 0, 'name', 'trucks');

% Add an objective function (with linear and quadratic terms)
total_penalty = gurobi.QuadExpr();
total_truck_cost = gurobi.LinExpr();

% Note the quadratic term!
for t = days
    total_penalty = total_penalty + penalty_coefficient * (demand(t{1}) - delivery(t{1}))^2;
end

for t = days
    total_truck_cost = total_truck_cost + truck_price * trucks(t{1});
end

MIQP.setObjective(total_penalty + total_truck_cost, 'Minimize');

% Add demand and capacity constraints
MIQP.addConstr(sum(cell2mat(values(delivery))) == sum(cell2mat(values(demand))));

for t = days
    MIQP.addConstr(delivery(t{1}) <= truck_size * trucks(t{1}));
end

% Optimize
MIQP.optimize();

% Display results
disp(['Total costs: ', num2str(MIQP.objVal)]);
disp(['Penalty costs: ', num2str(total_penalty.getValue())]);
disp(['Truck cost: ', num2str(total_truck_cost.getValue())]);

for t = days
    disp(['Deliver ', num2str(delivery(t{1}).X), ' (demand ', num2str(demand(t{1})), ') with ', num2str(trucks(t{1}).X), ' trucks']);
end