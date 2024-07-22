clc
clear
names = ['A','B','C','D','E']; %Same as in part c), but added again for completeness.
V_e = [3,2,1,7,6; 40,5,6,7,8; 4,9,5.5,10,6];
k = [2,3,1];

k_th_price(names, V_e, k)

function cell_array = k_th_price(names, values, k)
all_indices = zeros(length(names),length(k));
all_prices = zeros(1,length(k));
for i = 1:length(k)
    V = values(i,:);
    k_i = k(i);
    [index, price] = k_th_buyers(names,V,k_i);
    all_indices(:,i) = index;
    all_prices(1,i) = price;
end
end

function [indices, price] = k_th_buyers(names, values, k)
iter = 1;
indices = zeros(length(names),1);
while iter <= k
    [values_sorted,index] = sort(values,'descend');
    values(index(1)) = [0];
    indices(index(1),1) = 1;
    iter = iter + 1;
end
price = max(values);
end