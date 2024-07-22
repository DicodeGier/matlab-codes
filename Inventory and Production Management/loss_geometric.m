function y = loss_geometric(x, p)
% LOSS_GEOMETRIC Loss function for the geometric distribution
%
% y = loss_geometric(x, p)
%
% Computes the loss function for the geometric distribution, i.e.,
%
% y = E[max(X-d,0)]
%
% where X is geometrically distributed with parameter 'p'.

y = (1-p).^(x+1) ./ p;
