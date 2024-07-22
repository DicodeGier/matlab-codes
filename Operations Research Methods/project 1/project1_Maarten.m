clear
clc

C=[9999    3    5   48   48    8    8    5    5    3    3    0    3    5    8    8  5
    3 9999    3   48   48    8    8    5    5    0    0    3    0    3    8    8   5
    5    3 9999   72   72   48   48   24   24    3    3    5    3    0   48   48  24
   48   48   74 9999    0    6    6   12   12   48   48   48   48   74    6    6  12
   48   48   74    0 9999    6    6   12   12   48   48   48   48   74    6    6  12
    8    8   50    6    6 9999    0    8    8    8    8    8    8   50    0    0   8
    8    8   50    6    6    0 9999    8    8    8    8    8    8   50    0    0   8
    5    5   26   12   12    8    8 9999    0    5    5    5    5   26    8    8   0
    5    5   26   12   12    8    8    0 9999    5    5    5    5   26    8    8   0
    3    0    3   48   48    8    8    5    5 9999    0    3    0    3    8    8   5
    3    0    3   48   48    8    8    5    5    0 9999    3    0    3    8    8   5
    0    3    5   48   48    8    8    5    5    3    3 9999    3    5    8    8   5
    3    0    3   48   48    8    8    5    5    0    0    3 9999    3    8    8   5
    5    3    0   72   72   48   48   24   24    3    3    5    3 9999   48   48  24
    8    8   50    6    6    0    0    8    8    8    8    8    8   50 9999    0   8
    8    8   50    6    6    0    0    8    8    8    8    8    8   50    0 9999   8
    5    5   26   12   12    8    8    0    0    5    5    5    5   26    8    8 9999];

C = C(1:10,1:10);
k2 = 45;
k3 = 10;
E = 5;

%function [lambda, zld, zlr] = tspLD(C, K2, K3, E, w)
%note: output x is of the form [x11 x21 x31 ... ]
n = size(C,1);

ones_vector = ones(1,n);
first_constraint = kron(eye(n),ones_vector);
second_constraint = repmat(eye(n),1,n);
constraint_matrix = [first_constraint;second_constraint];

rhs = [ones(2*n,1);ones(k2,1);repmat(2,k3,1)];

preallocated_2 = zeros(k2, n*n);

all_combinations_2 = nchoosek(1:n,2);
for i = 1:k2
    indices = all_combinations_2(i,:);
    output_matrix = zeros(n);
    output_matrix(indices,indices) = 1;
    constraint_to_add = output_matrix(:);
    preallocated_2(i, :) = transpose(constraint_to_add);
end
constraint_matrix = [constraint_matrix; preallocated_2];

w = ones(n,1);
w(4,1) = 0;
w(7,1) = 0;

%Enforce S ⊆ W
loop_vec = 1:n;
loop_vec = loop_vec(w == 1);

preallocated_3 = zeros(k3, n*n);
subhandle = loopchoose(loop_vec,3);
to_add = zeros(k3, 3);
for i = 1:k3

    to_add(i,:) = subhandle();
    indices_3 = to_add(i,:);
    output_matrix_3 = zeros(n);
    output_matrix_3(indices_3,indices_3) = 1;
    constraint_to_add_3 = output_matrix_3(:);
    preallocated_3(i,:) = transpose(constraint_to_add_3);
end

constraint_matrix = [constraint_matrix; preallocated_3];

dualized_indices = zeros(E,size(loop_vec,2));
dualized_constraint_matrix = zeros(E,n^2);
dualized_rhs = zeros(E,1);
random_numbers = randi([4 size(loop_vec,2)-1],1,E);
random_numbers = sort(random_numbers);

for i = 1:E
    if i == 1 | random_numbers(i-1) ~= random_numbers(i)
        subhandle = loopchoose(loop_vec,random_numbers(1,i));
    end
    indices_4 = subhandle();
    dualized_indices(i,1:size(indices_4,2)) = indices_4;
    output_matrix_4 = zeros(n);
    output_matrix_4(indices_4,indices_4) = 1;
    constraint_to_add_4 = output_matrix_4(:);
    dualized_constraint_matrix(i,:) = transpose(constraint_to_add_4);
    dualized_rhs(i,1) = random_numbers(i)-1;
end

%Presets subgradient algorithm
%NOTE we take a negative lambda
lambda = -1 * ones(E , 1); 
stop = 0;
k = 1;
mu_zero = 1;
alpha = 0.9;
epsilon = 0.05;
%
%We can define the lagrangian relaxation outside of the loop with the
%objective function as the exeption as this is the only component changing
%with lambda
main_model.A = sparse(constraint_matrix);
main_model.rhs = rhs;
main_model.sense = [repmat('=',2*n,1);repmat('<',k2+k3,1)];
main_model.modelsense = 'min';
main_model.vtype = repmat('B',size(constraint_matrix,2),1);
params.outputflag = 0;
while stop == 0
    
    main_model.obj = transpose(C(:)) - transpose(lambda) * dualized_constraint_matrix;
    result = gurobi(main_model, params);
    
    xsol = reshape(result.x,[n n]);
    
    gamma_x = dualized_rhs - (dualized_constraint_matrix * result.x);
    
    if gamma_x == zeros(E,1)%1 element van lambda of alle elements??
        break
    end
    
    mu = mu_zero + (alpha^k);
    
    lambda_check = lambda + mu*gamma_x; %deze zit er voor nu ff bij om te kijken
    %bij development 
    lambda_next = min(lambda + mu*gamma_x,0);
    if norm(lambda_next - lambda) > epsilon
        k = k+1;
        lambda = lambda_next;
    else
        break
    end
    
end

%end


%to do for later
%prealocate constraint matrix
%rhs toevoegen met preallocation

%in de gaten houden
% k2 en k3 voor max loop iterations
%Geen output_matrix_3 nodig voor eindproduct