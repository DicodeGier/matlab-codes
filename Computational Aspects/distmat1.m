function d = distmat1(x)
%DISTMAT1 Distance matrix (one for-loop).
%
%   D = DISTMAT1(X) returns the distance matrix with all distances
%   between the points represented by the rows of X.
m = size(x, 1);
d = zeros(m, m);  % initialize output matrix
for i = 1:m-1
  % Take the difference of i-th row with all of the next rows.
  row_diff = x(i, :) - x(i+1:m, :);
  
  % Calculate the norm for all rows using elementwise arithmetic.
  dist = sqrt(sum(row_diff.^2, 2));
  
  % Store results (symmetric)
  d(i+1:m, i) = dist;
  d(i, i+1:m) = dist';
end
x = [2,3;4,5;6,7]
distmat1(x)