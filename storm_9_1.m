function stat_prob = storm_9_1(c, rho)
total = 0;
for n = 0:(c-1)
    total = total + ((c*rho)^n)/factorial(n);
end
stat_prob = 1/(1+(1-rho)*(factorial(c)/((c*rho)^c))*total);  
end

