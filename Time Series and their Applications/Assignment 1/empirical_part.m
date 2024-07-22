clear
clc

%load the data
data = load("data.mat");
co2 = data.co2;
ret = data.ret;

%%%Question a%%%
%mean
mean_ret = mean(ret);
%variance
var_ret = var(ret);
%ACF's
[ACF,lags] = autocorr(ret,20);

%%%Question b%%%
p_value = LB(ret,20);

%%%Question c%%%
%quantiles_plot = qqplot(ret);

%%%Question d%%%
abs_ret = abs(ret);
p_value_abs = LB(abs_ret,20);

%%%Question e%%%
n = length(co2);
R = 0;
for t = n:-1:1 % work backwards such that each t > s
    value_t = co2(t);
    for s = t-1:-1:1 % consider all values of s below t
    value_s = co2(s);
    if value_t > value_s
        R = R+1;
    end
    end
end
mu_R = n*(n-1)/4;
sigma_squared_R = n*(n-1)*(2*n+5)/72;

Z = (R - mu_R - 0)/(sqrt(sigma_squared_R)/sqrt(n));
 
lower_bound_CI = mu_R - (1.96*sqrt(sigma_squared_R)/sqrt(n));
upper_bound_CI = mu_R + (1.96*sqrt(sigma_squared_R)/sqrt(n));

%%%Question f%%%
beta = regress(co2, [ones(n,1) (1:n)']);

x = 1:n;
y =  beta(1) + beta(2)*x;
%{
plot(co2)
hold on
plot(x,y)
%}
residuals = co2' - y;
R_residual = 0;
for t = n:-1:1 % work backwards such that each t > s
    value_t_residual = residuals(t);
    for s = t-1:-1:1 % consider all values of s below t
    value_s_residual = residuals(s);
    if value_t_residual > value_s_residual
        R_residual = R_residual+1;
    end
    end
end

Z_residual = (R_residual - mu_R - 0)/(sqrt(sigma_squared_R)/sqrt(n));

function p = LB(x,h)
    T = length(x);
    [ACF_function, ~] = autocorr(x,h);
    ACF_function(1) = []; %remove rho_X(0) from the list
    total_sum = 0;
    for i = 1:h
        total_sum = total_sum + (ACF_function(i)^2)/(T-i);
    end
    Q_LB = T*(T+2)*total_sum;
    p = 1 - chi2cdf(Q_LB,h);        
end


