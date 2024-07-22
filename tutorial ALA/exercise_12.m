function S_m = exercise_12(m)
n = 1;
total = 2;
while n <= m
    product = (2*n)^2/((2*n)^2 - 1);
    total = total * product;
    n = n + 1;
end
disp(total)
end
