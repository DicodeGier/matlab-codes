                                                                                                                                                                                                                                                                                                        clear
clc
close all

data = load('data.mat');
data = data.sp500;
returns = data.Return;

sub_data = data(data.Date <= "2024-03-28",:);
sub_returns = sub_data.Return;
T = size(sub_returns,1);

%%
%%% 1 %%%
all_results = cell(36,6);

for p = 0:5
    for q = 0:5
        all_results{p*6+q+1,1} = p;
        all_results{p*6+q+1,2} = q;
        
        mdl = arima(p,0,q);
        result = estimate(mdl,sub_returns);
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

row_index = find([all_results{:,4}] == min([all_results{:,4}])); % determine 'best' based on AIC
opt_p = all_results{row_index,1};
opt_q = all_results{row_index,2};

opt_mdl = arima(opt_p,0,opt_q);
[opt_result,cov,logL,info] = estimate(opt_mdl,sub_returns);

%%
clc
%%% 2 %%%
res = infer(opt_result,sub_returns);
[~, pValue] = archtest(res);

res_squared = res.^2;
autocorr(res_squared)


%% 
clc
%%% 3 %%%
opt_mdl_garch = arima(opt_p,0,opt_q);
opt_mdl_garch.Variance = garch(2,2);
opt_result_garch = estimate(opt_mdl_garch,sub_returns);
summarize(opt_result_garch)

%%
clc
close all
%%% 4 %%%
h = size(data,1) - T;
garch_res_estimated = infer(opt_result_garch.Variance,sub_returns);
v = forecast(opt_result_garch.Variance,h,garch_res_estimated);

%%% 5 %%%
picture = plot(1:984, data.RV_Proxy(1:984));
ax = ancestor(picture, 'axes');
ax.YAxis.Exponent = 0;
hold on
plot(1:984,garch_res_estimated)
legend('Proxy', 'Estimated')
grid on
xlabel('Observations in time')
ylabel('Conditional variance')
title('Estimated Conditional Variance against Proxy')

hold off

picture_2 = scatter(1:22,v);
ax = ancestor(picture_2, 'axes');
ax.YAxis.Exponent = 0;
hold on 
scatter(1:22,data.RV_Proxy(985:end))
xlabel('s')
ylabel('Conditional Variance')
title('Prediction vs Proxy')
grid on
legend('Prediction', 'Proxy')
