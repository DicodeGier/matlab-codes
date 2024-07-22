function y = loss2_normal(x, mu, sigma)
% LOSS2_NORMAL Second-order loss function for the normal distribution
%
% y = loss2_normal(x, mu, sigma)
%
% Computes the second-order loss function for the normal distribution,
% i.e.,
%
% y = 0.5*E[max(X-x,0)^2]
%
% where X is normally distributed with mean 'mu' and standard deviation
% 'sigma'.

if nargin==1
  mu = 0;
  sigma = 1;
end

[x,mu,sigma] = samesize(x,mu,sigma);

z = (x-mu)./sigma;
y = sigma.^2 .* ((z.^2+1) .* (1-normcdf(z)) - z.*normpdf(z)) / 2;
