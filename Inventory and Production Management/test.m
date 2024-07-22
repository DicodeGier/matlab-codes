
L = [75 85 95]
test2 = [75 85 95 95]
array = []

for j = 1:size(L,2)
    pos = test2 == L(j)
    entry = sum(pos)
    array(:,end+1) = entry
end