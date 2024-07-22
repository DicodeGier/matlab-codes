clear

%% Initialize
N = 1000;                               % no. of replications
m = 1000;                               % no. of rows
n = 40;                                 % no. of columns

% Generate random matrices
X = rand(m, n);
A = rand(n, n);

%% Using a standard loop
% We replicate the experiment N times to get a more reliable timing
% comparison.
tic
for t = 1:N
    [m, n] = size(X);
    y = zeros(m, 1);
    for i=1:m
        y(i) = X(i,:) * A * X(i,:)';
    end
end
toc

%% Vectorized, but many unnecessary computations
tic
for t=1:N
    % Compute x(i,:) * A * x(j,:)' for all i and j (in a vectorized way)
    % and keep only those where i==j (diagonal of matrix).
    z = diag(X * A * X');
end
toc

%% Using an efficient vectorized approach
tic
for t=1:N
    % Note that the i-th row of X*A is X(i,:)*A. We multiply this result
    % elementwise by X(i,:). To get the desired result, we only need to
    % take the sum along the 2nd dimension.
    z = sum((X * A) .* X, 2);
end
toc
