function [errors, index_loc] = ex_6_6_4(func, der_func, x0, stepsize)
% ex_6_6_4 returns matrix with errors and minimum error location
%
% Input
% -----
%   func: function we want to evaluate
%   der_func: derivative of func
%   x0: point at which function is evaluated
%   stepsize: stepsize at horizontal axis
%
% Output
% ------
%   errors: matrix of errors
%   index_loc: index location at which minimum error is achieved
x0 = x0(:)';
stepsize = stepsize(:);

d = (func(x0 + stepsize) - func(x0 - stepsize)) ./ (2*stepsize);
errors = abs(der_func(x0) - d);
[~, index_loc] = min(errors);



end

