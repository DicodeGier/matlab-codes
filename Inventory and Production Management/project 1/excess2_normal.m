function y = excess2_normal(x, mu, sigma)
% EXCESS2_NORMAL Second-order excess function for the normal distribution
%
% y = excess2_normal(x, mu, sigma)
%
% Computes the second-order excess function for the normal distribution,
% i.e.,
%
% y = 0.5*E[max(x-Z,0)^2]
%
% where Z is normally distributed with mean 'mu' and standard deviation
% 'sigma'.

if nargin==1
  mu = 0;
  sigma = 1;
end

y = 0.5*(sigma.^2 + (x-mu).^2) - loss2_normal(x,mu,sigma);
