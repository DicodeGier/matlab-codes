clear

f = @(x) [1-x(1:2:end-1)', 10*(x(1:2:end-1)' - x(2:2:end)').^2];
rng default %Fix randomness of the script
n = 10;
x0 = rand(n,1);
option = optimoptions('fsolve','Display','iter-detailed');
sol = fsolve(f,x0,option);

T = 20
func_eval = @(T,x0) fsolve(f,x0,optimoptions('fsolve','Display','iter-detailed','MaxFunctionEvaluations',T));
sol2 = func_eval(T,x0)

