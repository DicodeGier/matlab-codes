function c = acf(phi, theta, h)
% This function computes the ACF of an ARMA process, 
% phi(L)X_t = theta(L)Z_t.
% The parameter h is the number of desired ACFs (r(1), ..., r(h)).
    
    c = acvf(phi, theta, h);
    c = c(2:end) / c(1);

end