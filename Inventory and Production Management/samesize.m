function varargout = samesize(varargin)
% SAMESIZE Converts all input arguments to the same size.
%
% [OUT1,OUT2,...] = SAMESIZE(IN1,IN2,...)
%
% Converts input arguments (scalar, vector or matrix) to output arguments
% which all have the same dimension.  Gives error if dimensions don't
% agree, except for dimension 1.

if nargin ~= nargout
  error('Number of input and output arguments must be equal.');
end
for n=1:nargin
  [r(n),c(n)] = size(varargin{n});
end
if any([r c] == 0)
  warning('samesize:emptydata', 'Empty input found.  All outputs are empty.');
  for i = 1:nargin
    varargout{n} = [];
  end
  return;
end
rows = max(r);
cols = max(c);

% Test for dimension 1 < dim < max_dim
r_test = all([r~=1; r~=rows],1);
c_test = all([c~=1; c~=cols],1);
if any([r_test c_test])
  error('Dimensions don''t agree');
end

% Dimensions are OK.
% Compute multiplication dimensions
r = (r==rows) + rows*(r==1 & r~=rows);
c = (c==cols) + cols*(c==1 & c~=cols);
for n=1:nargin
  varargout{n} = repmat(varargin{n},r(n),c(n));
end

