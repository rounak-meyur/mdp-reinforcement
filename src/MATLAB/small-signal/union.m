function [v] = union(v1,v2)
% Syntax: [v] = union(v1,v2)
%
% Purpose: Find the union of the integers in the vectors 
%          v1 and v2.
%
% Input: v1, v2 - vectors of integers
%
% Output: v - a vector containing all the integers in 
%             v1 and v2.
%
% See Also: V_slow, coh_map, L_group

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
if length(v1) >= length(v2)
    [nrow,ncol] = size(v1);
  else
    [nrow,ncol] = size(v2);
end
if nrow == 1
    vx = [v1 v2];
  else
    vx = [v1; v2];
end
vx = sort(vx);
nvx = length(vx)
v = vx(1); 
nv = 1;
for i = 2:n_vx 
  if vx(i) ~= v(nv)
    nv = nv + 1;
    v(nv) = vx(i);
  end
end
if nrow ~= 1
  v = v'
end
