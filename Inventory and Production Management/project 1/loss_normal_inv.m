function x = loss_normal_inv(y, mu, sigma)
% LOSS_NORMAL_INV Inverse normal loss function
%
% x = loss_normal_inv(z, mu, sigma)
%
% Returns the value of x for which
%
%   E[max(Z-x,0)] = y
%
% where Z is normal distributed with mean 'mu' and standard deviation
% 'sigma'.

if nargin==1
  mu = 0;
  sigma = 1;
end

[mu, sigma, y] = samesize(mu, sigma, y);

x = zeros(size(y));
for i=1:numel(y)
  f = @(z)(loss_normal(z,mu(i),sigma(i)) - y(i));
  x(i) = fzero(f, mu(i));
end
