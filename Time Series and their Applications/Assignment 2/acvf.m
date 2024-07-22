function c = acvf(phi, theta, h)
% This function computes the ACVF of an ARMA process, 
% phi(L)X_t = theta(L)Z_t with Var(Z_t) = 1.
% The parameter h is the number of desired ACVFs
% (c(1) is the variance of X).
    
    phi = phi(:); theta = theta(:);
    p = length(phi) - 1;
    q = length(theta) - 1;
    r = max(p, q);
    
    if r == 0
       c = [((theta / phi) ^ 2); zeros(h, 1)];
       return
    end
    
    h = h + 1;
    thi = flipud(theta);
    
    if q < r
       thi = [thi; zeros(r - q, 1)];
    end
    
    g = mpbf(thi, thi);
    
    if p == 0
       nm = min(q + 1 + h, length(g));
       c = g(q + 1 : nm);
       c = [c; zeros(h - length(c) + 1, 1)];
    else
       phii = phi;
       
       if p < r
          phii = [phii; zeros(r - p, 1)];
       end
       
       A = zeros(r + 1);
       
       for i = 1 : r + 1
           A(i, 1) = phii(i);
           A(i, r + 1) = phii(r - i + 2);
           
           if i > 1
              for k = 2 : i
                  A(i, k) = A(i, k) + phii(i - k + 1);
                  A(i, r + 2 - k) = A(i, r + 2 - k) + phii(r + 1 - i + k);
              end
           end
       end
    
       g = A \ g(1 : r + 1);
       c = polydiv(flipud(g), phii(1 : p + 1), h - 1)';
       c(1) = 2 * c(1);
    end

    c = c(1:h);

end

function x = mpbf(a, b)

    b = flipud(b);
    x = conv(a, b);

end
