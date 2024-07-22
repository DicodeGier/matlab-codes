a = 5;
b = 2;

f = @(x) gampdf(x,a,b);


expectation = integral(@(x) x.*f(x), 0, inf)
variance = integral(@(x) x.^2.*f(x), 0, inf) - expectation^2
covariance = integral2(@(x,y) x.*f(x).*y.*f(y), 0, inf, 0, inf) - expectation^2

