function x = sim_arch_m(mu, delta, alpha0, alphaq, T)

    b = 500;
    alphaq = flip(alphaq(:));
    p = numel(alphaq);
    z = zeros(T + b + p, 1);
    v = randn(T + b, 1);
    h = zeros(T + b, 1);
    x = zeros(T + b, 1);

    for t = 1:(T + b)
        h(t) = alpha0 + (z(t : t + p - 1) .^ 2)' * alphaq;
        z(t + p) = sqrt(h(t)) * v(t);
        x(t) = mu + delta * h(t) + z(t + p);
    end

    x = x((b + 1):end);
end