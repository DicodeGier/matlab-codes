clear
clc

names = ['A','B','C','D','E'];
V_c = [3,2,1,3,3; 4,5,6,7,8; 5,10,5,10,5];
V_c_cell = num2cell(V_c,2)
[indices, price] = cellfun(@assignment,[names,V_c_cell],'uniformoutput',false)

function [index, price] = second_price(v)
index = find(v == max(v), 1, 'last');
price = max(v([1:index-1, index+1:end]));
end

function overview = assignment(names, values)
[index, price] = second_price(values)
%%overview = struct('all_winners_per_item', all_winners, 'all_prices_per_item', all_prices);
end