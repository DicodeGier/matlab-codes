function y = loss_gamma(x, rho, lambda)
% LOSS_GAMMA Loss function for gamma distribution.
%
% y = loss_gamma(x, rho, lambda)
%
% Returns E[ max(Z-x, 0) ], where Z is gamma distributed with shape parameter
% 'rho' and scale parameter 1/'lambda'.

if nargin < 3
  lambda = 1;
end

[x,rho,lambda] = samesize(x,rho,lambda);

% Do some conversion so that we can call the gamma distribution with
% lambda=1, which is numerically more stable.
z = x.*lambda;
y = (rho.*(1-gamcdf(z,rho+1)) - z.*(1-gamcdf(z,rho))) ./ lambda;
