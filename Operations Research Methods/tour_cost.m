function cost = tour_cost(tour,A)
    %{
    %%%
    function calculates cost of a tour
    %%%
    input: tour of which you want to know the cost and distance matrix A
    %%%
    output: cost of the tour
    %}
    n = length(tour);
    cost = 0;
    for i=1:n-1
        cost = cost + A(tour(i),tour(i+1));
    end
    cost = cost + A(tour(n),tour(1));
end

