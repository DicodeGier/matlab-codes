clear
clc

mu_L = 230
sigma_L = 60
s = 350
Q = 600

%1
P2 = 1 - loss_normal(s,mu_L,sigma_L)/Q

%2
exp_x = 1*0.1+10*0.4+100*0.4+250*0.1
exp_x_squared = 1*0.1+100*0.4+10000*0.4+62500*0.1
exp_x_cubed = 1*0.1+1000*0.4+1000000*0.4+15625000*0.1

exp_u = exp_x_squared/(2*exp_x) - 0.5
sigma_u = sqrt(exp_x_cubed/(3*exp_x) - (exp_x_squared/(2*exp_x))^2 - 1/12)

%3
mu_w = mu_L + exp_u
sigma_w = sqrt(sigma_L^2 + sigma_u^2)

P2_undershoots = 1 - (loss_normal(s,mu_w,sigma_w)/Q)

%4
s_solver = @(s) ((loss_normal(s,mu_w,sigma_w) - loss_normal(s + Q,mu_L,sigma_L))/(1-P2_undershoots) - Q - exp_u)
s_small = fzero(s_solver,100)
S_capital = Q + s_small

%5
lambda = mu_w / sigma_w^2;
rho = mu_w^2 / sigma_w^2;
eqn = @(s)((loss_gamma(s, rho, lambda) - loss_gamma(s+Q, rho, lambda))/(1-P2_undershoots) - Q - exp_u);
s_gamma = fzero(eqn, mu_w);
S_gamma = s_gamma + Q
