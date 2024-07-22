clear
clc
%%%%%%%%%%%%%%%%%%%%%%%
%Input
n = 4;
W = 7;

c = [10,7,25,24];
w = [2,1,6,5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%checks for feasibility of the input
assert(n>=1, 'n should be at least 1')
assert(W>=0, 'W should be at least 0')
assert(all(c>0)&&all(w>0),'values and weights should be nonnegative')

%initialize output
output = zeros(W+1,n+1);
for m = 1:n
    if m ==1
        output(1:W+1-w(1),2) = c(1);
    else
    for W_prime = 1:W
        previous_column = output(W_prime,m);
        if W_prime <= W+1-w(m)
            adjustment = output(W_prime+w(m),m)+c(m);
            max_output = max(previous_column,adjustment);
            output(W_prime,m+1) = max_output;
        else
        output(W_prime, m+1) = previous_column;
        end
        
    end
    end  
end
full_matrix = output
optimal_value = output(1,n+1)




