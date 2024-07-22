%% Initialize
clear
m = 2000;         % number of points
n = 10;           % dimension
N = 10;           % no. of replications
X = randn(m, n);  % create random matrix

%% N replications for the solution with 2 loops
tic
for i=1:N
  D = distmat2(X);  % calculcate distance using 2 loops
end
toc

%% N replications for the solution with a single loop
tic
for i=1:N
  D = distmat1(X);  % calculcate distance using 1 loop
end
toc
