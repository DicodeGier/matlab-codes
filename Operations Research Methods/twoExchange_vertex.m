function [opt2_v,tour] = twoExchange_vertex(A,stour,sopt)
%This code searches along the 2-exchange vertex neighbourhood to find the
%optimal tour, starting from the tour 'stour'

%+-----------+-----------------------------------------------+
%| Input:    |                                               |
%+-----------+-----------------------------------------------+
%| A         | Weighted adjacency matrix of graph G          |
%| stour     | Initial tour in G                             |
%| sopt      | Cost of stour                                 |
%+-----------+-----------------------------------------------+
%| Output:   |                                               |
%+-----------+-----------------------------------------------+
%| opt2_v    | Cost of optimal 2-exchange vertex tour        |
%| tour      | Optimal 2-exchange vertex tour                |
%+-----------+-----------------------------------------------+

n = length(A); 
tour = stour; 
tmin = 0; 

a = tour(n); 
for i = 1:n-1
    a_pred = a; %a_pred is always the predecessor of vertex a in the tour
    a = tour(i);
    a_suc = tour(i+1); 
    for j = i+1:n
        b = tour(j);
        b_pred = tour(j-1); %b_pred is always the predecessor of vertex b in the tour
        b_suc = tour(mod(j,n)+1); %b_suc is always the successor of vertex b in the tour
        
        %To compute the cost change, we consider three cases:
        if a_suc == b 
            %b comes directly after a:
            test = A(a_pred, b) + A(a,b_suc) - A(a_pred,a) - A(b,b_suc);
        elseif b_suc == a
            %a comes directly after b:
            test = A(b_pred, a) + A(b,a_suc) - A(b_pred,b) - A(a,a_suc); 
        else
            %there is at least one vertex between a and b:
            test = A(a_pred, b) + A(b,a_suc) + A(b_pred,a) + A(a,b_suc) ...
                - A(a_pred,a) - A(a,a_suc) - A(b_pred,b) - A(b, b_suc);
        end
        
        
        %Save the best reduction:
        if test < tmin
            tmin = test; 
            imin = i; 
            jmin = j; 
        end
        
    end
end

%Update tour if tmin is smaller than 0. This means that the tour remains
%unchanged if tmin = 0.
if tmin < 0
    %Exchange the vertices on positions i and j in the tour:
    a = tour(imin);
    b = tour(jmin); 
    tour(imin) = b;
    tour(jmin) = a;
end

opt2_v = sopt + tmin; 

end