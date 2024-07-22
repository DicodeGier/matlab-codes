clear
clc
close all

data = load('data.mat');
data = data.data;
subsample = data(data.Date <= "2022-12-31",:);
T = size(subsample,1);

%%% 1 %%%
sub_temp = subsample.Temp;
plot(subsample.Date,sub_temp)
xlabel("Date")
ylabel("Temperature Deviation")
title("Temperature Deviation over time")
grid on
 
January = [];
February = [];
March = [];
April = [];
May = [];
June = [];
July = [];
August = [];
September = [];
October = [];
November = [];
December = [];
for i = 1:T
    date = subsample.Date(i);
    temperature = subsample.Temp(i);
    splitted = split(string(date),"-");
    month = splitted(2);
    if month == "01"
        January(end+1) = temperature;
    elseif month == "02"
        February(end+1) = temperature;
    elseif month == "03"
        March(end+1) = temperature;
    elseif month == "04"
        April(end+1) = temperature;
    elseif month == "05"
        May(end+1) = temperature;
    elseif month == "06"
        June(end+1) = temperature;
    elseif month == "07"
        July(end+1) = temperature;
    elseif month == "08"
        August(end+1) = temperature;
    elseif month == "09"
        September(end+1) = temperature;
    elseif month == "10"
        October(end+1) = temperature;
    elseif month == "11"
        November(end+1) = temperature;
    elseif month == "12"
        December(end+1) = temperature;
    end
    
end

q25 = [quantile(January,0.25) quantile(February,0.25) quantile(March,0.25) quantile(April,0.25) quantile(May,0.25) quantile(June,0.25) quantile(July,0.25) quantile(August,0.25) quantile(September,0.25) quantile(October,0.25) quantile(November,0.25) quantile(December,0.25)];
q50 = [quantile(January,0.50) quantile(February,0.50) quantile(March,0.50) quantile(April,0.50) quantile(May,0.50) quantile(June,0.50) quantile(July,0.50) quantile(August,0.50) quantile(September,0.50) quantile(October,0.50) quantile(November,0.50) quantile(December,0.50)];
q75 = [quantile(January,0.75) quantile(February,0.75) quantile(March,0.75) quantile(April,0.75) quantile(May,0.75) quantile(June,0.75) quantile(July,0.75) quantile(August,0.75) quantile(September,0.75) quantile(October,0.75) quantile(November,0.75) quantile(December,0.75)];

plot(1:12,q25)
hold on
plot(1:12,q50)
hold on
plot(1:12,q75)
xlabel('Month')
ylabel("Temperature Deviation")
title("Quantiles for Temperature Deviation per Month")
grid on
legend("25th Quantile", "Mean", "75th Quantile")
hold off

%%% 2 %%%
D12 = LagOp({1 -1},'Lags',[0,12]);
trans_sub_temp = filter(D12,sub_temp);

autocorr(trans_sub_temp)
parcorr(trans_sub_temp)

[h,pValue] = adftest(trans_sub_temp);

%%% 4 %%%
test_data = data(877:end,:);
load("data_modeler.mat")
periods = size(data,1) - T;
[forecasts, MSE] = forecast(SARIMA_sub_temp15,periods,sub_temp);

scatter(test_data.Date, test_data.Temp)
hold on
scatter(test_data.Date,forecasts)
xlabel('h')
ylabel('Value')
title('Forecasts against Realized Values')
legend('Realized', 'Forecasted')
grid on
hold off

%%%%%%%%%%%%%%%%%%%%%%%
%%%     part 2      %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% 1 %%%
TSI = data.TSI;
DFT = fft(TSI);
DFT_conj = conj(DFT);
periodogram = DFT.*DFT_conj;

N = size(periodogram,1);
frequencies = (0:N-1) / N;
idx = frequencies > 0 & frequencies <= 0.5;
selected_frequencies = frequencies(idx);
selected_periodogram = periodogram(idx);

plot(selected_frequencies,selected_periodogram)
xlabel('Fourier Frequencies')
ylabel('Values')
title('Periodogram of TSI against Fourier Frequencies')
grid on

plot(selected_frequencies,log(selected_periodogram))
xlabel('Fourier Frequencies')
ylabel('Log Values')
title('Periodogram of Log TSI against Fourier Frequencies')
grid on

%%% 2 %%%
Weights = [3/33,3/33,3/33,3/33,3/33,2/33,1/33];
f_hat_0 = Weights(1)*periodogram(2) + 2* Weights(2:7)*periodogram(3:8);

output = [];
omegas_all = 0:0.001:0.5;
for omegas = omegas_all
    f_hat = DSAE(omegas,Weights,frequencies,periodogram,f_hat_0);
    output(end+1) = f_hat;
end

plot(selected_frequencies,selected_periodogram)
hold on
plot(omegas_all,output)
xlabel('Fourier Frequencies')
ylabel('Values')
title('Raw Periodogram vs Smoothed Periodogram')
legend('Raw Periodogram', 'Smoothed Periodogram')
grid on
hold off

plot(selected_frequencies,log(selected_periodogram))
hold on
plot(omegas_all,log(output))
xlabel('Fourier Frequencies')
ylabel('Log Values')
title('Log Raw Periodogram vs Log Smoothed Periodogram')
legend('Log raw Periodogram', 'Log smoothed Periodogram')
grid on
hold off

%%% 3 %%%
mdl_AR = arima('ARLags',[1 2 3 5 6 9 11 18 21]);
mdl_estimated = estimate(mdl_AR,TSI);
phi = cell2mat(mdl_estimated.AR);
theta = 1;
sigma2 = mdl_estimated.Variance;

implied_spd = arma_spd([1 -phi],theta,sigma2,omegas_all);

plot(omegas_all,output)
hold on
plot(omegas_all,implied_spd)
xlabel('Frequencies')
ylabel('Value')
title('Implied Spectral Density vs Discrete Spectral Average Estimate')
legend('Discrete spectral average estimate', 'Implied spectral density')

hold off

plot(omegas_all,log(output))
hold on
plot(omegas_all,log(implied_spd))
xlabel('Frequencies')
ylabel('Log Value')
title('Log Implied Spectral Density vs Log Discrete Spectral Average Estimate')
legend(' Log discrete spectral average estimate', 'Log implied spectral density')
grid on

function f_hat = DSAE(omega,Weights,frequencies,periodogram,f_hat_0)
    
    periodogram(1) = f_hat_0;
    periodogram = [flip(periodogram(end-5:end));periodogram];
    f_hat = 0;
    for i = -6:6
        f_hat = f_hat + Weights(abs(i)+1)*periodogram((find(min(abs(frequencies-omega))==abs(frequencies-omega),1))+i+6);    
    end
end