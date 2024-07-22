function [y,c] = ex_3_11_5(x)
    d = diff(x);
    i = [find(d), length(x)];        % always include the last entry
    y = x(i);
    c = zeros(1,length(d));
    counter = 0;
    index = 1;
    for j = 1:length(x)
        if j == 1
            value = x(1);
            counter = counter + 1;
        else
            if x(j) == value && j == length(x)
                c(index) = counter + 1;
            elseif x(j) ~= value && j == length(x)
                c(index) = counter
                index = index + 1
                c(index) = 1
            elseif x(j) == value
                counter = counter + 1;  
            else
                value = x(j);
                c(index) = counter;
                counter = 1;
                index = index + 1;
            end
        end
    end
    k = c>0;
    c = c(k)
end