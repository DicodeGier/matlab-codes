clear
n = 30
x = 4
S(1) = x^0/factorial(0) + x^1/factorial(1)
k = 2
a = (x^k)/(factorial(k))
iter = 0

while abs(a)>0.0001 & iter<= n-1
    iter = iter + 1
    a = (x^k)/(factorial(k))
    S(k) = S(k-1) + a
    k = k + 1
    if iter == n-1
        disp('more than 30 terms are added')
    end
end
