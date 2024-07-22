clc
clear all;
close all;

%% Group details
% Group 21 consisting of
% 1. Daan Spanbroek
% 2. Daan van Turnhout
% 3. Dico de Gier
% 4. Hendrik (Nick) Verkaik
%% Excercise 1c (run function of part b) on the input below)
names = ['A','B','C','D','E'];
V_c = [3,2,1,3,3; 4,5,6,7,8; 5,10,5,10,5];
struct_array = assignment(names, V_c)
%% Exercise 1e (run function of part d) on the input below)
names = ['A','B','C','D','E']; %Same as in part c), but added again for completeness.
V_e = [3,2,1,7,6; 40,5,6,7,8; 4,9,5.5,10,6];
k = [2,3,1];
cell_array = k_th_price(names, V_e, k)
%% Exercise 1h) (run function of part g) on the input below)
G = {makedist('Gamma','a',5), makedist('Gamma','a',4), ...
        makedist('Gamma','a',3)};
n = 10;
r = [3,6,4];
expected_revenue = total_expected_revenue(G,r,n)
%% Exercise 1i) (experimental verification)
G_2 = {makedist('halfnormal','mu',0,'sigma',2)};
n = 50;
k = 1;
expected_revenue_halfnormal = total_expected_revenue(G_2,k,n)

pd = G_2{1}; %take the above object for the random function
rng(1) %%set seed to control for randomness
v = random(pd,100000,n);
number_of_rows = size(v,1);
total_price = 0;
for i = 1:size(v,1)
    [~,price] = second_price(v(i,:));
    total_price = total_price + price;
end
average_price_charged = total_price/number_of_rows

%% Exercise 1k) (running function of part j) on some input)
values = gamrnd(5,1,100000,2);
p_tilde = 1:10;

output = second_price_reserve(values,p_tilde); %calculate the values
mean_output = mean(output); %take the mean
plot(p_tilde,mean_output,'-o','MarkerIndices',1:length(p_tilde))
title('average profit against reserve price')
xlabel('reserve price')
ylabel('average profit')
%%now determine the maximum point
ymax = max(mean_output);
ymax_index = find(mean_output == ymax);
xmax = p_tilde(ymax_index);
hold on
plot(xmax,ymax,'square','MarkerSize',13);
hold off
legend('all points','maximum point')

%% Function of Exercise 1a)
function [index, price] = second_price(v)
index = find(v == max(v), 1, 'last'); %find the index of maximum price
price = max(v([1:index-1, index+1:end])); %find maximum price without the above index
end
%% Function of Exercise 1b)
function S = assignment(names,V)
%generates a structure array, runs SECOND_PRICE on every row of V 
    S = struct('Item',{},'Name',{},'price',{}); %initialize struct
    for n = 1:size(V,1)
        [i,b] = second_price(V(n,:)); %for every row of V, apply function of a) 
        S(n).Item = n; S(n).Name = names(i); S(n).price = b; %store the results in the struct
    end
end

%% Function of Exercise 1d)
function S = k_th_price(names, V, k)
    %initialize arrays to store the results
    receivers = zeros(size(V));
    prices = zeros(1,length(k));
    for n = 1:length(k)
        [~,i] = maxk(V(n,:),k(n)); %returns indices of the k highest bids in the row of V
        receivers(n,i) = 1; %store which buyers get a copy
        %set bids of the above buyers to zero and determine the max ...
        %value (= price) of the remaining values in the row of V
        V(n,i) = 0;
        price = max(V(n,:));
        prices(n) = price;
    end
S = cell(length(names),3);
%now put the results in a cell_array
for j = 1:length(names)
    S(j,1) = {names(j)};
    S(j,2) = {receivers(:,j)};
    S(j,3) = {prices*receivers(:,j)}; %vector product to determine total price
end
end
%% Function of Exercise 1f)
function anon_func  = function_handle_creater(object,k,n)
anon_func = @(x)(factorial(n)./(factorial(k).*factorial(n-k-1)).*pdf(object,x).*(cdf(object,x)).^(n-k-1).*(1-cdf(object,x)).^k);
end
%% Function of Exercise 1g)
function total_revenue = total_expected_revenue(cell_array,k,n)
total_revenue = 0;
for i = 1:length(k)
    %take the corresponding elements
    object = cell_array{i};
    k_i = k(i);
    %create the function handle
    output = function_handle_creater(object,k_i,n);
    %create new function to calculate the expectation
    new_func = @(y) (y.*output(y));
    %integrate over all values to find the expectation
    expectation = integral(new_func,-Inf,Inf);
    %see also derivations in report as to why the following step is valid
    revenue = expectation * k_i;
    total_revenue = total_revenue + revenue;
end
end
 
%% Function of Exercise 1j)
function price = second_price_reserve(values,p_tilde)
dim = size(values,1);
price = zeros(dim,length(p_tilde));
for i = 1:dim
    for j = 1:length(p_tilde)
    row_values = values(i,:);%take the row of the matrix
    highest_price = max(row_values);%check what the highest bid is
        if highest_price < p_tilde(j) %if highest bid < reserve price --> price = 0
            price(i,j) = 0;
        else
          %determine the second highest price
          [~,index] = sort(row_values,'descend');
          row_values(index(1)) = [0];
          second_max = max(row_values);
          %price is max of reserve price and second highest price
          price(i,j) = max(p_tilde(j), second_max);
        end
    end
end
end