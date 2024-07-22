clear
n = [4, 2, 10]
A = 10 * rand(n)
i = (A>5) & (A<6)
A(i) = round(A(i))
A(~i) = round(2*A(~i))/2
A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = find(A>=6)
B = zeros(length(ind),3)
[B(:,1), B(:,2), B(:,3)] = ind2sub(n, ind) %%vragen, documentatie ind2sub zei niks over 3-dim
B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
average_grade = mean(reshape(A,[8,10]))