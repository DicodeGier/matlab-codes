function y = loss2_poisson(x, mu)
% LOSS2_POISSON Second order loss function for the Poisson distribution
%
% y = loss2_poisson(x, mu)
%
% Computes the second-order loss function for the poisson distribution,
% i.e.,
%
% y = 0.5 * E[max(X-x,0)*max(X-x-1,0)]
%
% where X is Poisson(mu) distributed.

y = (mu.*(mu-x).*poisspdf(x,mu) + ((mu-x).^2+x).*(1-poisscdf(x,mu)))/2;
