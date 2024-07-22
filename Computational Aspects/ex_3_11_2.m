clear
A = rand(1000,10);
max_row = size(A,1);
distance_matrix_1 = zeros(max_row,max_row);
tic
for i = 1:max_row
    for j = 1:max_row
        if i==j
            continue
        else
            distance_matrix_1(i,j) = sqrt(sum((A(i,:) - A(j,:)).^2));
        end
    end
end
two_loop = toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = rand(1000,10);
max_row = size(A,1);
distance_matrix_2 = zeros(max_row,max_row);
tic
for i = 1:max_row
  difference = A(i,:) - A(i+1:max_row,:);
  entry = sqrt(sum(difference.^2, 2));
  distance_matrix_2(i,i+1:max_row) = entry;
  distance_matrix_2(i+1:max_row,i) = entry;
end
one_loop = toc



    