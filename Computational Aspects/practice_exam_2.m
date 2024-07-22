clear
clc

fraction_larger = @(x) (sum(x>mean(x))/length(x));
x = [9, 19, 16, 20, 14, 1, 17, 19];
f = fraction_larger(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 20;
n = 5;
X = rand(m, n);
y = rand(m, 1);

beta_ols = X \ y
%%%
Aeq = ones(1,n);
beq = 1;
beta_second_method = lsqlin(X,y,[],[],Aeq,beq,zeros(1,n),[]);
%%%
SSE_1 = (y - X*beta_ols)'*(y-X*beta_ols)
SSE_1_check = sum((y-X*beta_ols).^2)
SSE_2 = (y - X*beta_second_method)'*(y-X*beta_second_method)
SSE_2_check = sum((y-X*beta_second_method).^2)
perc_increase = (SSE_2 - SSE_1)/SSE_1 * 100
perc_increase_check = (SSE_2_check - SSE_1_check)/SSE_1_check * 100
