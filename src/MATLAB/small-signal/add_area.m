function [area,narea] = add_area(area,narea,i,k)
% Syntax: [area,narea] = add_area(area,narea,i,k)
%
% Purpose: Add machine i to area k
%
% Input: area - a matrix of machines in the areas.
%        narea - a vector of the number of machines in each 
%                area.
%        i - machine number.
%        k - area number.
%
% Output: area - new area
%         narea - new number of machines in each area
%
% See Also: add_area

% Algorithm: 
%
% Calls:
%
% Call By: 

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
n = length(narea);
if k <= n  % add to existing area
  temp = sort([area(k,1:narea(k)) i]);
  area(k,1:narea(k)+1) = temp;
  narea(k) = narea(k) + 1;
 else      % create new area
  area(k,1) = i;
  narea(k,1) = 1;
end  

