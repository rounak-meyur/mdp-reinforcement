function [k] = element(x,v)
% Syntax: [k] = element(x,v)
%
% Purpose: Test if x is an element of the vector v
%
% Input: x - integer
%        v - a vector of integers
%
% Output: k - 0 if x is not an element of v
%             1 otherwise.
%
% See Also: 

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
n_v = length(v); 
k = 0;
diff = min(abs(x*ones(size(v))-v));
if round(diff) == 0
  k = 1;
end
