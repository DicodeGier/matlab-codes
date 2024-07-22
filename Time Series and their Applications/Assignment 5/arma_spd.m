function f = arma_spd(phi, theta, sigma2, v)
% This function computes the spectral density of an ARMA process, 
% phi(L)X_t = theta(L)Z_t, with Var(Z_t) = sigma2

    T = polyval(flip(theta), exp(-2 * pi * 1i * v));
    F = polyval(flip(phi), exp(-2 * pi * 1i * v));
    f = sigma2 * abs(T) .^ 2 ./ abs(F) .^ 2;

end