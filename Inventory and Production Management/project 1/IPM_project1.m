%% data 
clear
clc

format longg

%load data
data = readtable("dataset06.csv")
data = table2array(data)

%%%extra parameters given in introduction
%fixed ordering cost:
A = 50
%r unknown, but proportional to v
%if order is out of stock, customers have to wait until new replenishment arrives
%backorder cost per euro value short
B2 = 0.50
%forecasting done weekly
%demand over different periods is i.i.d.

%%%extract parameter from table for later use
id = data(:,1)
value = data(:,2)
expected_demand = data(:,3)
stdev_demand = data(:,4)
leadtime = data(:,5)
Q = data(:,6)

%replenish item when current inventory position reaches projected demand during lead time + 1 extra week as buffer
s = expected_demand.*leadtime + expected_demand

%% Replenishment policy
%check for assumption of normal distribution: if v = sigma/mu>0.5 use Gamma
normal_check = stdev_demand./expected_demand
normal_check_binary = normal_check > 0.5
number_normal = 25 - sum(normal_check_binary)
number_gamma = sum(normal_check_binary)
fprintf('# items with normal distribution is %g and #items with gamma distribution is %g',[number_normal, number_gamma])
mean_v = mean(normal_check)
fprintf('mean coefficient of variation is %.4f which is larger than 0.5 => gamma should be used\n',[mean_v])

%check for ignorance of shortages at start of RC, start again with normal
mu_l = expected_demand.*leadtime
sigma_l = stdev_demand.*sqrt(leadtime)
k = (s-mu_l)./sigma_l
ESPRC = loss_normal(s,mu_l,sigma_l)
ESPRC_corrected = loss_normal(s,mu_l,sigma_l) - sigma_l.*loss_normal(k+Q./sigma_l)

difference = abs(ESPRC-ESPRC_corrected)
mean_difference = mean(difference)

%check for undershoots
%done in words

%% Current performance
rho = (mu_l.^2)./(sigma_l.^2);
lambda = (sigma_l.^2)./mu_l;
%average on hand stock
on_hand_stock = Q./2 + s - mu_l
%p1 service level
p1 = round(gamcdf(s,rho,lambda),4)
%p2 service level
%we want 1/lambda for loss_gamma
lambda = lambda.^(-1);
p2 = round(1-(loss_gamma(s,rho,lambda) - loss_gamma(s+Q,rho,lambda))./Q,4)


%% EOQ
N = sum(expected_demand./Q)
TACS = sum(0.5.*Q.*value)
r = Q.^-2 .*2 .* A .*expected_demand ./ value

exchange = @(N) (0.5 * (sum(sqrt(expected_demand.*value)))^2) ./ N;fplot(exchange)

hold on
plot([0,N],[TACS,TACS],"k--")
hold on
plot([N,N],[0,TACS],"k--")

grid on
xlim([5,25])
ylim([100000,200000])
title("EOQ-exchange curve")
xlabel("N")
ylabel("TACS")
hold on
plot(N,TACS,"r.",'MarkerSize',20)
text(N+0.5,TACS+1400,"Current policy")

N_opt = 0.5 * (sum(sqrt(expected_demand.*value)))^2 / TACS;
value_2 = TACS/N_opt;
hold on
plot(N_opt,TACS,"r.",'MarkerSize',20)
text(N_opt-3.2,TACS-2400,"EOQ policy")

r = N_opt/TACS*A
order_quantities = sqrt(2.*A.*expected_demand./(value.*r))
order_quantities_floored = floor(order_quantities./10)*10
order_quantities_ceiled = ceil(order_quantities./10)*10

TRC_lower = A.*expected_demand./order_quantities_floored + 0.5.*order_quantities_floored.*value.*r
TRC_upper = A.*expected_demand./order_quantities_ceiled + 0.5.*order_quantities_ceiled.*value.*r
lower_better = TRC_lower < TRC_upper
final_order_quantities = order_quantities_floored.*lower_better + order_quantities_ceiled.*(~lower_better)

%% fill rate
p2_target = 0.98;
lambda = (sigma_l.^2)./mu_l;
%we want 1/lambda for loss_gamma
lambda = lambda.^(-1);
rho = (mu_l.^2)./(sigma_l.^2);

%p2 service equation
eqn = @(s) loss_gamma(s,rho,lambda) - loss_gamma(s+final_order_quantities,rho,lambda) - final_order_quantities.*(1-p2_target);
s_new = fsolve(eqn,mu_l);

s_new = ceil(s_new)
% check if p2 >= 0.98
new_p2 = round(1-(loss_gamma(s_new,rho,lambda) - loss_gamma(s_new+final_order_quantities,rho,lambda))./final_order_quantities,4)

%% eval cost
lambda = (sigma_l.^2)./mu_l;
rho = (mu_l.^2)./(sigma_l.^2);

order_cost_current = 52.*A.*expected_demand./Q
order_cost_fillrate = 52.*A.*expected_demand./final_order_quantities

on_hand_stock_current = Q./2 + s - mu_l;
on_hand_stock_fillrate = final_order_quantities./2 + s_new - mu_l;
holding_cost_current = 52.*on_hand_stock_current.*value.*r
holding_cost_fillrate = 52.*on_hand_stock_fillrate.*value.*r

k_fillrate = (s_new-mu_l)./sigma_l;
ESPRC_gamma = @(s,Q) loss_gamma(s,rho,lambda.^(-1)) - loss_gamma(s+Q,rho,lambda.^(-1));
shortage_cost_current = 52.*B2.*value.*expected_demand./Q .*ESPRC_gamma(s,order_quantities)
shortage_cost_fillrate = 52.*B2.*value.*expected_demand./final_order_quantities .*ESPRC_gamma(s_new,final_order_quantities)

total_cost_current = order_cost_current + holding_cost_current + shortage_cost_current
total_cost_fillrate = order_cost_fillrate + holding_cost_fillrate + shortage_cost_fillrate

%% implied cost
lambda = (sigma_l.^2)./mu_l;
B2_new = (final_order_quantities.*r)./(expected_demand.*(1-gamcdf(k_fillrate,rho,lambda)))
k_check = gaminv(1-((final_order_quantities.*r)./(expected_demand.*B2_new)),rho,lambda)
s_check = k_check.*sigma_l + mu_l
p2_check = 1 - (loss_gamma(s_check,rho,lambda.^(-1)) - loss_gamma(s_check + final_order_quantities,rho,lambda.^(-1)))./final_order_quantities

%% minimize cost
%we apply the B2 decision rule
lambda = (sigma_l.^2)./mu_l
eqn_s = @(s) value.*r + B2.*value.*expected_demand./final_order_quantities.*(-rho./lambda.*gampdf(s,rho+1,lambda) + gamcdf(s,rho,lambda) + s.*gampdf(s,rho,lambda) + rho./lambda.*gampdf(s+final_order_quantities,rho+1,lambda) - gamcdf(s+final_order_quantities,rho,lambda) - (s+final_order_quantities).*gampdf(s+final_order_quantities,rho,lambda))
s_min = fsolve(eqn_s,s)

order_cost_min = A.*expected_demand./final_order_quantities
holding_cost_min = (final_order_quantities./2+s_min-mu_l).*value.*r
shortage_cost_min = B2.*value.*ESPRC_gamma(s_min,final_order_quantities).*expected_demand./final_order_quantities
total_cost_min = order_cost_min + holding_cost_min + shortage_cost_min
 
