c = addme(1,2)

function b = addme(a1,a2,a3)
    if nargin == 1
        b = a1;
    elseif nargin == 2
        b = a1 + a2;
    else
        b = a1 + a2 + a3;
    end
end