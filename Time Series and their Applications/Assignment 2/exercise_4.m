clear

%%% exercise 4 %%%
Gamma = full(gallery('tridiag', 4, -0.5, 1.25, -0.5));
gamma = [-0.5;0;0;0];
MSPE_vector = [];
a_t_cell = {};

for i = 1:4
    a_t = Gamma(1:i,1:i) \ gamma(1:i);
    a_t_cell{1,i} = a_t; 

    MSPE = 1.25 - a_t'*gamma(1:i);
    MSPE_vector = [MSPE_vector;MSPE];
end