function [B] = diagb(A,area,n_area)
% Syntax: [B] = diagb(A,area,n_area)
%
% Purpose: Set the off-diagonal blocks of the matrix A equal %          to 0.
%
% Input: A - a square matrix
%        area - a matrix defining the block diagonal 
%               structure of A.
%        n_area - a vector containing the number of induces
%                 in the diagonal blocks.
%
% Output: B - A with diagonal blocks equal to 0.
%
% See Also: offdb

% Algorithm: 
%
% Calls: offdb
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

B = offdb(A,area,n_area);
B = A - B;
