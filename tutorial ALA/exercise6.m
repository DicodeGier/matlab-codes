clear
C = [2,4,4,20; 1,5,8,22; 0,2,4,8]
C(1,:) = C(1,:)/2
C(2,:) = C(2,:)- C(1,:)
C(2,:) = C(2,:)/3
C(1,:) = C(1,:) - 2*C(2,:)
C(3,:) = C(3,:) - 2*C(2,:)
D = rref(C)

%the third column does not have a pivot, hence there is one free variable.
%Therefore there are infinitely many solutions in the form [2+2*x3; 4-2*x3;
%x3]

% no since there is a free variable, there are infinitely many solutions to
% the system Ax = 0 but linear independence requires that only the trivial
% solution is a solution to Ax = 0

p = [2;4;0]
v = [2;-2;1]

A = [2,4,4; 1,5,8; 0,2,4]
b = [20; 22; 8]
x = A\b