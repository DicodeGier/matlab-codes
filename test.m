clear
clc
%test
p = [2 4 2 0 6 4 6];
treeplot(p)
[x,y] = treelayout(p);
text(x+0.02,y,{1,2,3,4,5,6,7})

                

    