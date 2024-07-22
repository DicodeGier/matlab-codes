            weight1 = parent_weights{1,j};
            weight2 =  parent_weights{2,i};
            childweight = (weight1 + weight2)/2;
            %{
            if min_matrices == max_matrices
                childweight = (weight1 + weight2)/2;
            else
                childweight = round(W_max/max_matrices) - 1;
            end
            %}
            if mod(childweight,1) ~= 0
                childweight = round(childweight) - 1;
            end