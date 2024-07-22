x = [1,3;2,4];
A = [1,2;3,4];

product = x.' * A * x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 1000;   
n = 40;
repetitions = 1000;

X = rand(m,n);
A = rand(n,n);

tic
for counter = 1:repetitions
    [m,n] = size(X);
    y = zeros(m,1);
    for i=1:m
      y(i) = X(i,:) * A * X(i,:)';
    end
end
toc

tic
for counter = 1:repetitions
    product = X * A * X';
end
toc
