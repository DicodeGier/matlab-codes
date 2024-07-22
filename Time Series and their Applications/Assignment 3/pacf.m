function p = pacf(phi, theta, h)
% This function computes the PACF of an ARMA process, 
% phi(L)X_t = theta(L)Z_t.
% The parameter h is the number of desired PACFs (a(1), ..., a(h)).
    
    c = acvf(phi, theta, h);
    n = length(c) - 1;
    f = zeros(n, 1);
    f(1) = c(2) / c(1);
    p = f;
    v = (1 - f(1) ^ 2) * c(1);
    
    if n > 1
        for i = 2 : n
            f(i) = (c(i + 1) - f(1:(i - 1))' * flipud(c(2:i))) / v;
            p(i) = f(i);
            v = (1 - f(i) ^ 2) * v;
            f(1:(i - 1)) = f(1:(i - 1)) - f(i) * flipud(f(1:(i - 1)));
        end
    end

end