function y = loss_normal(x, mu, sigma)
% LOSS_NORMAL Loss function for normal distribution
%
% y = loss_normal(x, mu, sigma)
%
% Returns E[ max(Z-x, 0) ], where Z is normally distributed with mean 'mu'
% and standard deviation 'sigma'.

if nargin==1
  mu = 0;
  sigma = 1;
end

[x,mu,sigma] = samesize(x,mu,sigma);

k = (x-mu)./sigma;
y = sigma .* (normpdf(k) - k .* normcdf(-k));
