%{
preallocated_2 = zeros(k2, n*n);
subhandle = loopchoose(1:n,2);
to_add_2 = zeros(k2,2);
for i = 1:k2
    to_add_2(i,:) = subhandle();
    indices = to_add_2(i,:);
    output_matrix = zeros(n);
    output_matrix(indices,indices) = 1;
    constraint_to_add = output_matrix(:);
    preallocated_2(i, :) = transpose(constraint_to_add);
end
constraint_matrix = [constraint_matrix; preallocated_2];

preallocated_3 = zeros(k3, n*n);
subhandle = loopchoose(loop_vec,3);
to_add_3 = zeros(k3, 3);
for i = 1:k3

    to_add_3(i,:) = subhandle();
    indices_3 = to_add_3(i,:);
    output_matrix_3 = zeros(n);
    output_matrix_3(indices_3,indices_3) = 1;
    constraint_to_add_3 = output_matrix_3(:);
    preallocated_3(i,:) = transpose(constraint_to_add_3);
end

constraint_matrix = [constraint_matrix; preallocated_3];
%}

sizes = randi([4 size(loop_vec,2)-1],1,E);
sizes= sort(sizes);