function y = loss_poisson(x,mu)
% LOSS_POISSON - First-order loss function for the Poisson distribution
%
% y = loss_poisson(x,mu)
%
% Computes the loss function for the poisson distribution, i.e.,
%
% y = E[max(X-y,0)]
%
% where X is Poisson(mu) distributed.

y = (mu-x).*(1-poisscdf(x,mu)) + mu.*poisspdf(x,mu);
