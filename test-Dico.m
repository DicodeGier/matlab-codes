a = 1:5
b = [2,1,3,4,5]
c = sort([2,1,3,4,5])

d = [1,2,3;4,5,6;7,8,9]
r = 2
e = d([1:r-1,r+1:end],[1:r-1,r+1:end])
minimum = min(min(e))
[row,col] = find(e == minimum)