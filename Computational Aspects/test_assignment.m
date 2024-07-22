clc
clear

values = gamrnd(5,1,10000,2);
p_tilde = 1:10;

output = second_price_reserve(values,p_tilde);
mean_output = mean(output);
plot(p_tilde,mean_output,'-o','MarkerIndices',1:length(p_tilde))
title('average profit against reserve price')
xlabel('reserve price')
ylabel('average profit')
ymax = max(mean_output);
ymax_index = find(mean_output == ymax);
xmax = p_tilde(ymax_index);
hold on
plot(xmax,ymax,'square','MarkerSize',13);

function price = second_price_reserve(values,p_tilde)
dim = size(values,1);
price = zeros(dim,length(p_tilde));
for i = 1:dim
    for j = 1:length(p_tilde)
    row_values = values(i,:);
    highest_price = max(row_values);
        if highest_price < p_tilde(j)
            price(i,j) = 0;
        else
          [~,index] = sort(row_values,'descend');
          row_values(index(1)) = [0];
          second_max = max(row_values);
          price(i,j) = max(p_tilde(j), second_max);
        end
    end
end
end

