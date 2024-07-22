clear
clc

all_K = [3 3 2];
number_routes = 3;

binary_matrix = zeros(number_routes, sum(all_K));
row = 1;
column = 1;

for n = 1:length(all_K)
    element = all_K(n);
    counter = 1;
    for operations = 1:element
        if counter <= element
            binary_matrix(row, column) = 1;
            column = column + 1;
            counter = counter + 1;
        end
    end
    row = row + 1;
end
 