function y = loss2_geometric(x, p)
% LOSS2_GEOMETRIC Second-order loss function for the geometric distribution
%
% y = loss2_geometric(x, p)
%
% Computes the second-order loss function for the geometric distribution,
% i.e.,
%
% y = E[ max(X-d,0) * max(X-d-1,0)]
%
% where X is geometrically distributed with parameter 'p'.

y = (1-p).^(x+2) ./ p.^2;
