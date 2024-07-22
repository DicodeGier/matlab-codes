clear
clc
close all
%%
data = load("data.mat");
data = data.data;

%%% (1) %%%
sub = data(data.Date <= "2022-08-31",:);
subspread = sub.Spread;
T = size(subspread,1);
all_results = cell(36,6);

for p = 0:5
    for q = 0:5
        all_results{p*6+q+1,1} = p;
        all_results{p*6+q+1,2} = q;
        
        mdl = arima(p,0,q);
        result = estimate(mdl,subspread);
        var = result.Variance;
        all_results{p*6+q+1,3} = var;
        
        log_var = log(var);
        dim = (1+p+q);
        AIC = log_var + dim*2/T;
        BIC = log_var + dim*log(T)/T;
        HQC = log_var + dim*2*log(log(T))/T;
        
        all_results{p*6+q+1,4} = AIC;
        all_results{p*6+q+1,5} = BIC;
        all_results{p*6+q+1,6} = HQC;
    end
end
%%
row_index = find([all_results{:,4}] == min([all_results{:,4}])); % determine 'best' based on AIC
opt_p = all_results{row_index,1};
opt_q = all_results{row_index,2};

opt_mdl = arima(opt_p,0,opt_q);
[opt_result,cov,logL,info] = estimate(opt_mdl,subspread);

%%
%%% 2 %%%

phi = cell2mat(opt_result.AR);
theta = cell2mat(opt_result.MA);

phi_input = [1 -phi];
theta_input = [1 theta];

theo_ACF = acf(phi_input,theta_input,20);
sample_ACF = autocorr(subspread);
sample_ACF = sample_ACF(2:end);
comparision_ACF = [theo_ACF sample_ACF];
theo_PACF = pacf(phi_input,theta_input,20);
sample_PACF = parcorr(subspread);
sample_PACF = sample_PACF(2:end);
comparison_PACF = [theo_PACF sample_PACF];

close all
scatter(0:20,[1;theo_ACF])
hold on
scatter(0:20,[1;sample_ACF])
grid on
xlabel('Lag')
ylabel('Autocorrelation')
title('Theoretical vs Sample ACF')
legend('Theoretical ACF', 'Sample ACF')

close all
scatter(0:20,[1;theo_PACF])
hold on
scatter(0:20,[1;sample_PACF])
grid on
xlabel('Lag')
ylabel('Autocorrelation')
title('Theoretical vs Sample PACF')
legend('Theoretical PACF', 'Sample PACF')

%%
%%% 3 %%%

summarize(opt_result)

res = infer(opt_result,subspread);

m = 20;
total_sum = 0;
sample_ACF_LB = autocorr(res,m);
sample_ACF_LB = sample_ACF_LB(2:end);
for i = 1:m
    total_sum = total_sum + (sample_ACF_LB(i)^2)/(T-i);
end
Q_LB = T*(T+2)*total_sum;
p = 1 - chi2cdf(Q_LB,m-dim);

%%
%%% 4 %%%
[forecasts, MSPE] = forecast(opt_result,7,subspread);
comparision_forecasts = [forecasts data.Spread(273:end)];

%%% 5 %%%
cb_upper = forecasts + 1.96 * sqrt(MSPE/T);
cb_lower = forecasts - 1.96 * sqrt(MSPE/T);

close all
picture = scatter(1:7,forecasts,'filled','blue');
ax = ancestor(picture, 'axes');
ax.YAxis.Exponent = 0;
grid on
hold on
scatter(1:7,data.Spread(273:end),'filled','green')
hold on
plot(1:7,cb_lower, 'red')
hold on
plot(1:7,cb_upper,'red')
xlabel('h')
ylabel('value')
title('Forecasts, Realizations and Confidence Bands')
legend('Forecasts','Realizations', 'Confidence Bands','')

