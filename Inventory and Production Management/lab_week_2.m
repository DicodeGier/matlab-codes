clear
clc

data = readtable('Lab1-2.csv');
data = table2array(data);

%%1
N = sum(data(:,1)./data(:,3));
TACS = sum(data(:,3)./2.*data(:,2));

%%2
exchange = @(N) (0.5 * (sum(sqrt(data(:,1).*data(:,2))))^2) ./ N;

%%3
fplot(exchange)
grid on
xlim([0,2000])
ylim([0,500000])
title("EOQ-exchange curve")
xlabel("N")
ylabel("TACS")
hold on;
plot(N,TACS,"r.",'MarkerSize',20)
text(N + 40, TACS + 10000, "Current point")

%%4
%%point beneath the current point
hold on;
value_1 = exchange(N)/N;
plot(N,exchange(N),"b.",'MarkerSize',20)
text(N + 40, exchange(N) + 10000, "A/r = "+round(value_1,2))
hold on
plot([0,N],[TACS,TACS],"--k")
hold on
plot([N,N],[0,TACS],"--k")

N_opt = 0.5 * (sum(sqrt(data(:,1).*data(:,2))))^2 / TACS;
value_2 = TACS/N_opt;
hold on
plot(N_opt,TACS,"b.",'MarkerSize',20)
text(N_opt - 30, TACS + 20000, "A/r = "+round(value_2,2))


%%5
red_N = (N_opt - N)/N * 100;
red_TACS = (exchange(N) - TACS)/TACS * 100;
check = red_N == red_TACS;



