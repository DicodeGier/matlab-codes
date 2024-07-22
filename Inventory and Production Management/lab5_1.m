%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is the template file for exercise 4.1 %%
%%                                            %%
%% Fill the missing parts marked by XXX       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Exercise SPT 9.3: Initialization

% Acquisition cost
v = 60;

% Selling price
p = 100;

% Salvage value
g = 51;

% Penalty cost beyond lost sales
B = 0;

% Discrete demand distribution
x = [300; 400; 500; 600; 700; 800];
px = [0.1; 0.1; 0.4; 0.2; 0.1; 0.1];


%% Moments of demand

% Derive mean and stdev from `x` and `px`.
mu = sum(x .* px);
sigma = sqrt(sum((x-mu).^2 .* px));

fprintf('Mean demand: %g\n', mu);
fprintf('Stdev demand: %g\n', sigma);


%% Order quantity discrete demand model

% Underage/overage cost
c_u = p - v + B;
c_o = v - g;

% Target newsvendor probability
target = c_u/(c_u + c_o);

% Cumulative probability distribution
Fx = cumsum(px);

% Smallest Q such that Fx >= target. Use find(., 1) to get this.
i = find(Fx >= target,1);                                % index smallest Q
Q = x(i);                               % value smallest Q

fprintf('\nTarget probability: %.4f\n', target);
fprintf('%6s %6s\n', 'Q', 'F(Q)');
fprintf('%6g %6.1f\n', [x, Fx]');
fprintf('Optimal Q = %g\n', Q);


%% Expected profit

% Anonymous function for profit for given Q and D
profit = @(Q, D) ((p-v)*D - c_o*max(Q-D,0) - c_u*max(D-Q,0));

% Anonymous function for expected profit (use `profit` from above)
exp_profit = @(Q) (sum(px .* profit(Q, x)));

% Initialize profit vector.
profit_Q = zeros(size(x));

% Compute expected profit for all possible order quantities
fprintf('\n');
fprintf('%6s %10s\n', 'Q', 'E(profit)');
for k = 1:length(x)
    profit_Q(k) = exp_profit(x(k));
    fprintf('%6d %10g\n', [x(k), profit_Q(k)]);
end
fprintf('\n');

% Plot order quantity vs. profit.
figure(1);
clf;
plot(x, profit_Q, '.-', 'MarkerSize', 20, 'Linewidth', 2);
xlabel('Order quantity');
ylabel('Expected profit');


%% Using normal distribution

% Newsvendor solution for normal distribution (unrounded).
Qn = norminv(c_u/(c_u+c_o),mu,sigma);

% Round to multiples of 100.
Qnr = round(Qn / 100) * 100;

% Print results
fprintf('Optimal order quantity for normal distribution:\n');
fprintf('Q = %.2f\n', Qn);
fprintf('Rounded to multiples of 100: %d\n\n', Qnr);


%% Exercise SPT 9.4: Discount

% Acquisition cost 55 if Q >= 750
v = 55;

% Updated underage/overage cost
c_u = p-v+B;
c_o = v-g;

% Target newsvendor probability
target = c_u/(c_u+c_o);

% Smallest Q such that Fx >= target (use find(., 1) to get this).
i = find(Fx >= target,1);                                % index smallest Q
Q = x(i);                               % value smallest Q

fprintf('Acquisition cost 55 if Q >= 750\n');
fprintf('\nTarget probability: %.4f\n', target);
fprintf('%6s %6s\n', 'Q', 'F(Q)');
fprintf('%6g %6.1f\n', [x, Fx]');
fprintf('Optimal Q = %g\n', Q);
if Q >= 750
    fprintf(['Unconstrained order point meets discount requirement,\n', ...
             'so Q=%d is now optimal.\n'], Q);
else
    fprintf(['Warning: optimal order quantity does not meet discount ' ...
             'requirement.\n']);
    % Now we have to compare the expected profit for Q=800 *with discount* with
    % the expected profit for Q=700 *without discount*.

end
