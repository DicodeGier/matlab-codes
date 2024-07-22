function NN_tour = initial_tour(A)
    n = size(A,1);
    NN_tour = [1];
    B = A;
    B(:,1) = 10000000;
    while length(NN_tour) < n
        minimum = min(B(NN_tour(end),:));
        new_vertex = find(B(NN_tour(end),:) == minimum);
        new_vertex = new_vertex(1);
        column_to_delete = new_vertex;
        NN_tour(end+1) = new_vertex;
        B(:,column_to_delete) = 10000000;
    end
end

