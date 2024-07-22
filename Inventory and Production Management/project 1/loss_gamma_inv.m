function x = loss_gamma_inv(y, rho, lambda)
% LOSS_GAMMA_INV Inverse gamma loss function
%
% x = loss_gamma_inv(y, rho, lambda)
%
% Returns the value of x for which
%
%   E[max(Z-x,0)] = y
%
% where Z is gamma distributed with shape parameter 'rho' and scale parameter
% 1/'lambda'.

if nargin < 3
  lambda = 1;
end

[y,rho,lambda] = samesize(y,rho,lambda);
x = zeros(size(y));
x0 = rho./lambda;
for i=1:numel(y)
  f = @(s)(loss_gamma(s,rho(i),lambda(i)) - y(i));
  x(i) = fzero(f, x0(i));
end
