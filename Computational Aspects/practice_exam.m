clear
clc
%%a
h = @(x,q,r,a,b) (gampdf(x,q,r)./(gamcdf(b,q,r) - gamcdf(a,q,r)));

%%b
g = @(x,q,r,a,b) (h(x,q,r,a,b).*(a<=x&x<=b));

%%c
x = 3:0.1:5;
plot(x,g(x,1,1,3,5))
xlabel('x');
ylabel('g');
title('values of g for x between 3 and 5');
hold on
plot(x,g(x,4,2,3,5))
hold on
plot(x,g(x,3,20,3,5))
legend('q = 1, r = 1', 'q = 4, r = 2', 'q = 3, r = 20')


%%
clear
clc

v = [1,5,2,3,4,9];
n = 3;
result = vector_placement(v,n)

function output_matrix = vector_placement(v,n)
    output_matrix = zeros(n,n);
    for i = 1:n
        input_vector = zeros(1,n);
        input_vector(:,end-i+1:end) = v(:,end-i+1:end);
        v = v(:,1:end-i);
        output_matrix(:,end-i+1) = input_vector';
    end     
end



