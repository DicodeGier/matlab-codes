clear
C = [2,3,1,21; 0,2,1,12; 1,1,2,18]

C(3,:) = C(3,:) - 0.5*C(1,:)
C(3,:) = C(3,:) + 0.25*C(2,:)

%the system has a unique solution since the rightmost column of the
%augmented matrix is not a pivot column, the solution is also unique, since
%there is no free variable

C(1,:) = C(1,:)/2
C(2,:) = C(2,:)/2
C(3,:) = C(3,:)/1.75
C(1,:) = C(1,:) - 1.5*C(2,:)
C(1,:) = C(1,:) + 0.25*C(3,:)
C(2,:) = C(2,:) -0.5*C(3,:)
D = rref(C)

x = [3; 3; 6]

A = [2,3,1; 0,2,1; 1,1,2]
b = A*x
x = A\b