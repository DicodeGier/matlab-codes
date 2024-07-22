function A = combsum_2(n,k)
x = cell(1,k);
[x{:}] = ndgrid(0:n);

A = zeros(numel(x{1}),k); %%numel counts number of elements
for i = 1:k
    A(:,i) = x{i}(:);
end
output = sum(A,2);
output = output == n;
A = A(output,:)
end

