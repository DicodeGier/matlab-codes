function length_MST = MST(W,r)
    number_vertices = size(W,1);
    W_helper = W; %copy W for later purposes
    T = [r]; %vector of vertices that are already added to the tree
    length_MST = 0; 
    W_outgoing = W(T,:); %W_outgoing will contain all feasible edges
    counter = 1;
    while counter < number_vertices
        if counter == 1
           %determine minimum length and add it to total length
           min_length = min(min(W_outgoing));
           length_MST = length_MST + min_length;
           %find index and add to T
           [row,col] = find(W == min_length);
           T(end+1) = setdiff(col, T);
           %update W, W_outgoing and counter
           W(:,col) = [];
           W_outgoing(:,col) = [];
           counter = counter + 1;
        else
        %first update W_outgoing with row of point that was added the
        %latest
        W_outgoing(end+1,:) = W(T(end),:);
        %determine min length, find index and update T
        min_length = min(min(W_outgoing));
        length_MST = length_MST + min_length;
        [row,col] = find(W_helper == min_length);
        T(end+1) = setdiff(col, T);
        %find the columns in W and W_outgoing that have to be deleted
        column_tobe_deleted = W_helper(:,T(end));
        for i = 1:size(W,2)
            if W(:,i) == column_tobe_deleted
                W(:,i) = [];
                W_outgoing(:,i) = [];
                break
            end
        end
        counter = counter + 1;
        end
    end
end
