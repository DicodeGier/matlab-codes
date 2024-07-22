possibilities = [];

for i = 0.1:0.1:1
    for j = 0.1:0.1:1
        for k = 0.1:0.1:1
            for l = 0.1:0.1:1
                for m = 0.1:0.1:1
                    for n = 0.1:0.1:1
                        if j + l + n == 1 && i + l + m == 1 && k + j + m == 1
                           possibilities(end+1,:) = [i,j,k,l,m,n];
                        end
                    end
                end
            end
        end
    end    
end