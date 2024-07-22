clear
clc

matrix = [Inf 3 4 4; 3 Inf 2 5; 4 2 0 6; 4 5 6 0];
source = 1;
destination = 4;
number_routes = 12;

[route, cost] = kShortestPath(matrix, source, destination, number_routes);