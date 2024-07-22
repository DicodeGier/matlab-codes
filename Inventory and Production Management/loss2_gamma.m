function y = loss2_gamma(x, rho, lambda)
% LOSS2_GAMMA Second-order loss function for the gamma distribution
%
% y = loss2_gamma(x, rho, lambda)
%
% Computes the second-order loss function for the gamma distribution,
% i.e.,
%
% y = 0.5*E[max(Z-x,0)^2]
%
% where Z is gamma distributed with shape parameter 'rho' and scale
% parameter '1/lambda'.

if nargin < 3
  lambda = 1;
end

[x,rho,lambda] = samesize(x,rho,lambda);

z = lambda.*x;
y = +rho.*(rho+1)./lambda.^2.*(1-gamcdf(z, rho+2))...
    -2*x.*rho./lambda.*(1-gamcdf(z, rho+1))...
    +x.^2.*(1-gamcdf(z,rho));
y = y/2;
