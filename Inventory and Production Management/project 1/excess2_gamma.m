function y = excess2_gamma(x, rho, lambda)
% EXCESS2_GAMMA Second-order excess function for the gamma distribution
%
% y = excess2_gamma(x, rho, lambda)
%
% Computes the second-order excess function for the gamma distribution,
% i.e.,
%
% y = 0.5*E[max(x-Z,0)^2]
%
% where Z is gamma distributed with shape parameter 'rho' and scale
% parameter '1/lambda'.

y = (rho.*(rho+1)./lambda.^2 - 2*x.*rho./lambda + x.^2)/2 - ...
    loss2_gamma(x, rho, lambda);

% Note that in general:
%
% excess2(d) = ( var(x) + (mu-d)^2 ) / 2 - loss2(d)
