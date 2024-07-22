clear
func = @sin;
der_func = @cos;
x0 = [0.5, 1, 2, -0.2];
len_x0 = length(x0);
n = 13;
h = 10.^-(0:n-1)';

[errors, index_loc] = ex_6_6_4(func, der_func, x0, h);

best_errors = zeros(len_x0,1);
for i = 1:len_x0
    best_errors(i,1) = errors(index_loc(i),i);
end

figure;
loglog(h, errors);
hold on
scatter(h(index_loc), best_errors, 'fill');
xlabel('step size');
ylabel('error');
title('Plot of Error vs. Step Size');
grid


