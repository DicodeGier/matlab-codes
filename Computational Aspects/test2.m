clc
clear

names = ['A','B','C','D','E']; %Same as in part c), but added again for completeness.
V_e = [3,2,1,7,6; 40,5,6,7,8; 4,9,5.5,10,6];
k = [2,3,1];
cell_array = k_th_price(names, V_e, k)

function S = k_th_price(names, V, k)
    receivers = zeros(size(V));
    prices = zeros(1,length(k));
    for n = 1:length(k)
        [~,i] = maxk(V(n,:),k(n));
        receivers(n,i) = 1;
        V(n,i) = 0;
        price = max(V(n,:));
        prices(n) = price;
    end
S = cell(length(names),3);
for j = 1:length(names)
    S(j,1) = {names(j)};
    S(j,2) = {receivers(:,j)};
    S(j,3) = {prices*receivers(:,j)};
end
end