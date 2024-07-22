clear
A = rand(3,10)
a1 = A>0.25 & A<0.5
a2 = all(A>0.25 & A<0.5)
a3 = A<0.25 | A>0.75
a4 = any(A<0.25 | A>0.75)