% ********************************************************************
%
% This is the template file for exercise 3_1
%
% Fill the missing parts marked by XXX
%
% ********************************************************************

clear

%% Input parameters
L  = 2.1;
R  = 1;
P2 = 0.96;


%% Normal distribution
mu_R    = 9.90;
sigma_R = 7.06;

% Derive parameters normal distribution for D(L) and D(R+L)
mu_L     = L/R*mu_R
sigma_L  = sqrt(L/R*(sigma_R^2))
mu_RL    = (1+L/R)*mu_R
sigma_RL = sqrt((1+L/R)*(sigma_R^2))

% Compute ESPRC assuming a normal distribution as a function of S.
% Include a correction for the shortage at the start of the RC.
ESPRC_n = @(S) loss_normal(S,mu_RL,sigma_RL) - loss_normal(S,mu_L,sigma_L);

% Define P2 service equation we need to solve as a function of S
eqn_normal = @(S) ESPRC_n(S) - (1-P2) * mu_R;

% Solve service equation using mu_RL as an initial guess
S_n = round(fzero(eqn_normal, mu_RL));
fprintf('Order-up-to level (normal case): S = %g\n', S_n);


%% Gamma distribution
mu_R    = 9.64;
sigma_R = 7.19;

% Derive parameters gamma distribution
rho_R  = (mu_R/sigma_R)^2
lambda = mu_R/(sigma_R)^2
rho_L  = rho_R*L/R
rho_RL = rho_R*(1+L/R)

% Compute ESPRC assuming a gamma distribution as a function of S
ESPRC_g = @(S) loss_gamma(S,rho_RL,lambda) - loss_gamma(S,rho_L,lambda);

% Define P2 service equation
eqn_gamma = @(S) ESPRC_g(S) - (1-P2) * mu_R;

% Solve service equation
S_g = round(fzero(eqn_gamma, mu_RL));
fprintf('Order-up-to level (gamma case): S = %g\n', S_g);


%% True service level
mu_R    = 11;
sigma_R = sqrt(50);

% Parameters true gamma distribution
rho_R  = (mu_R/sigma_R)^2
lambda = mu_R/(sigma_R)^2
rho_L  = rho_R*L/R
rho_RL = rho_R*(1+L/R)

% Redefine the ESPRC for the gamma distribution using the true parameters.
ESPRC_true = @(S) loss_gamma(S,rho_RL,lambda) - loss_gamma(S,rho_L,lambda);

% Calculate the true service levels (with the gamme distribution) for the two
% order-up-to levels
P_n = 1 - ESPRC_true(S_n)/mu_R
P_g = 1 - ESPRC_true(S_g)/mu_R

fprintf('True P2 service level for S=%g: %.4f\n', S_n, P_n);
fprintf('True P2 service level for S=%g: %.4f\n', S_g, P_g);
