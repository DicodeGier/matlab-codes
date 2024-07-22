function y = polydiv(b, a, n)
% psi(z) = psi0 + psi1 * z + ... = (b0 + b1 * z + ...) / (a0 + a1 * z + ...)
    
    x = [1 zeros(1, n)];
    y = filter(b, a, x);

end