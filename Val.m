function V = Val(k,b)
% N - number of different items
% w - array of weights for each item
% r - array of values for each item
m = floor(b/w(k)); % determine max number of item k for budget b
p = 0:m; % array of possible numbers of each item given budget b
if k==N
      V = max(r(k)*p); % base case
else
      temp = zeros(1,length(p));
      % recursion step
      for n=1:length(p)
          % value of k+1 item given budget remaining if p number of item k is
          % used
          temp(n) = Val(k+1,b-w(k)*p(n)); 
      end
      V = max(r(k)*p + temp); 
end
end