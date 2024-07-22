%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the key file for exercise 3.2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input data

D = [1000; 500; 100];                   % Yearly demand
v = [1.0; 0.8; 2.0];                    % Item value
Q = [200; 200; 100];                    % Order quantity
sigma = [100; 100; 30];                 % StDev demand during LT
TSSV = 120;                             % total safety stock value

% A quick reference to the local function (see below) capturing all fixed
% parameters.
my_summary = @(k) print_summary(k, D, v, Q, sigma);

% Compute the total safety stock value for a given safety factor.
tssv = @(k) sum(k .* sigma .* v);


%% Same P1

% Same P1 implies same safety factor k:
k = TSSV / sum(sigma .* v);
k = repmat(k, length(D), 1);
header('Same P1');
my_summary(k);


%% Same P2

% Derive the required safety factor for a given P2 target
safety_factor = @(P2) loss_normal_inv(Q .* (1-P2) ./ sigma);

% Search the P2 target that meets the TSSV target and compute corresponding
% safety targets.
eqn = @(P2) tssv(safety_factor(P2)) - TSSV;
P2 = fzero(eqn, 0.95);
k = safety_factor(P2);
header('Same P2');
my_summary(k);


%% Same B1

% Derive the required safety factor for a given ratio of B1 over r. By
% taking the maximum with 1 for the argument of the log we ensure that the
% safety is zero when the argument is less than 1.
arg = @(B1_over_r) B1_over_r * D / sqrt(2*pi) ./ (Q .* v .* sigma);

safety_factor = @(B1_over_r) sqrt(2*log(max(arg(B1_over_r), 1)));

% Now find the right ratio of B1 over r just as we did for the same P2.
eqn = @(B1_over_r) tssv(safety_factor(B1_over_r)) - TSSV;
B1_over_r = fzero(eqn, 100);
k = safety_factor(B1_over_r);
header('Same B1');
fprintf('Use B1/r = %g\n\n', B1_over_r);
my_summary(k);


%% Same B2

% Derive the required safety factor for a given ratio of B2 over r. If the
% argument is less than 0.5, then we force the safety factor to be equal
% to zero.
arg =  @(B2_over_r) 1 - Q ./ D / B2_over_r;
safety_factor = @(B2_over_r) norminv(max(0.5, arg(B2_over_r)));

% Now find the right ratio of B2 over r just as we did for the same P2.
eqn = @(B2_over_r) tssv(safety_factor(B2_over_r)) - TSSV;
B2_over_r = fzero(eqn, 1);
k = safety_factor(B2_over_r);
header('Same B2');
fprintf('Use B2/r = %g\n\n', B2_over_r);
my_summary(k);


%% Equal time supply

safety_factor = @(t) t * D ./ sigma;

% Now find the time supply t just as we did for the same P2.
eqn = @(t) tssv(safety_factor(t)) - TSSV;
t = fzero(eqn, 1);
k = safety_factor(t);
header('Equal time supply');
fprintf('Safety stock in units of time: %g\n\n', t);
my_summary(k);


%% Local functions in script file
function print_summary(k, D, v, Q, sigma)
% PRINT_SUMMARY summarize results for given safety factors.
%
% Here we derive the basic performance measures for all items as a function
% of the safety factors used.

SS = k .* sigma;
SSV = SS .* v;
P1 = normcdf(k);
ESPRC = sigma .* loss_normal(k);
P2 = 1 - ESPRC ./ Q;
ETSOPY = (1-P1) .* D ./ Q;
ETVSPY = ESPRC .* D ./ Q .* v;

fprintf([repmat(' %9s', 1, 8), '\n'], ...
        'k', 'SS', 'SSV', 'P1', 'P2', 'ESPRC', 'ETSOPY', 'ETSOPY');
fprintf(repmat('-', 1, 80));
fprintf('\n');
fprintf([repmat(' %9.4f', 1, 8), '\n'], ...
        [k, SS, SSV, P1, P2, ESPRC, ETSOPY, ETVSPY]');
fprintf(repmat('-', 1, 80));
fprintf('\n');
fprintf([repmat(' ', 1, 20), ' %9.4f', repmat(' ', 1, 30), ' %9.4f %9.4f\n']', ...
    sum(SSV), sum(ETSOPY), sum(ETVSPY));
fprintf('\n\n');

end


function header(title)
% HEADER Print a header line with title.

width = 80;
n = length(title);
line = repmat('*', 1, width);
k = floor((width - n - 4) / 2);
padder = repmat('*', 1, k);
fprintf('%s\n', line);
fprintf('%s%s%s\n', padder, pad(title, width - 2 * k, 'both'), padder);
fprintf('%s\n', line);

end
