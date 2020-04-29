function [A] = offdb(A,area,n_area)
% Syntax: [A] = offdb(A,area,n_area)
%
% Purpose: Set the diagonal blocks of the matrix A equal to 
%          0.
%
% Input: A - a square matrix
%        area - a matrix defining the block diagonal 
%               structure of A.
%        n_area - number of indices in diagonal blocks.
%
% Output: A - the input A with diagonal blocks equal to 0.
%
% See Also: diagb

% Algorithm: 
%
% Calls: 
%
% Call By: ex_group

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     May 1991
n = length(n_area); 
for i = 1:n
  ii = n_area(i);
  v = area(i,1:ii)';
  for j = 1:ii
    A(v,area(i,j)) = zeros(ii,1);
  end
end
