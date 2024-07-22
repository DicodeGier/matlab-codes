clear
clc
close all

data = load("data.mat");
data = data.data;

%%% (1) %%%
sub = data(data.Date <= "2022-08-31",:);

plot(sub.Date,sub.Spread)
xlabel("Date")
ylabel("Spread")
title("Monthly Australian interest rate spread")
legend("Spread")

autocorr(sub.Spread)
[ACF, lags] = autocorr(sub.Spread);

%%% (2) %%%
spread = sub.Spread;
mean_spread = mean(spread);
demeaned_spread = spread - mean_spread;

[result,~,~,~,statistics] = regress(demeaned_spread, lagmatrix(demeaned_spread, 1:4));

%%% (3) %%%
theoretical_acvf = acvf([1;-result],1,20);
theoretical_acf = theoretical_acvf./theoretical_acvf(1);
comparison = [theoretical_acf ACF];

plot(0:20,theoretical_acf)
hold on
plot(0:20,ACF)
xlabel('h')
ylabel('ACF')
title('Comparison of Theoretical and Sample ACF')
legend('Theoretical ACF', 'Sample ACF')
grid on

hold off

%%% (4) %%%
spread_2 = spread;
for h = 1:7    
    n = length(spread_2);
    forecast = sum(result.*spread_2(n:-1:n-3));
    spread_2(end+1,:) = forecast;
end

forecasts = spread_2(273:end);
realizations = data.Spread(273:end);
difference = realizations - forecasts;
var_wn = statistics(4);

MSPEs = [];
for h = 1:7
    MSPEs(h) = var_wn*sum((polydiv([1],[1 -result'], h-1).^2));
end

%%% (5) %%%
predictions = forecasts + mean_spread;
confidence_band = 1.96*sqrt(MSPEs/size(demeaned_spread,1));

picture = scatter(1:7,predictions,"filled");
ax = ancestor(picture, 'axes');
ax.YAxis.Exponent = 0;
hold on
plot(1:7,predictions' + confidence_band,"red")
hold on
plot(1:7,predictions' - confidence_band,"red")
hold on
scatter(1:7,realizations,"filled")
xlabel('h')
ylabel('value')
title('Forecasts with Confidence Bands and Realizations')
legend('Predictions', 'Confidence Bands','' ,'Realizations')
grid on

