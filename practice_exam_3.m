greaterthan = @(x,z) (sum(x>z));
x = [2,3,4,5,6];
z = 3;
greaterthan(x,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
mean_pos = @(x) (mean(x(x>0)));
x = [-1,0,1,2,3];
mean_pos(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
checkerboard([1,2,3,4,5;6,7,8,9,10;11,12,13,14,15])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
checkerboard_2([1,2,3,4,5,6;7,8,9,10,11,12;13,14,15,16,17,18;19,20,21,22,23,24])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
x = [1,2,3];
loop1(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
x = [1,2,3;4,5,6];
i = any(x<3,2)
x = x(~i,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load practice_exam
S = struct('word',{},'total',{});
words = unique(data(:,2));
for i = 1:length(words)
    word = words{i};
    j = strcmp(data(:,2),word);
    S(i).word = word;
    S(i).total = sum(cell2mat(data(j,1)));
end

function y = loop1(x)
n = length(x);
y = zeros(n-1,1);

i = 1:n-1;
y(i) = (2.*(x(i)+x(i+1))).*(x(i)+x(i+1)>=0) + -3.*(x(i)+x(i+1)).*(x(i)+x(i+1)<0);
end


function x = checkerboard_2(x)

[m,n] = size(x);
i = 1:2:m;
j = 1:2:n;
x(i,j) = 0;

i_2 = 2:2:m;
j_2 = 2:2:n;
x(i_2,j_2) = 0
end

function x = checkerboard(x)

[m,n] = size(x);
if mod(m,2) == 0
    x = 'error, number of rows is not odd'
else
    x(1:2:end) = 0
end
end