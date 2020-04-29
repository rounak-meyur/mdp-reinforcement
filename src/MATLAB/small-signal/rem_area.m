function [area,narea] = rem_area(area,narea,i)
% Syntax: [area,narea] = rem_area(area,narea,i)
%
% Purpose: Remove machine i from its area
%
% Input: area - a matrix of machines in the areas.
%        narea - a vector of the number of machines in each %                area.
%        i - machine number.
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

k = 0;
j = 0;
while k == 0
  j = j+1;
  temp = area(j,1:narea(j));
  k = element(i,temp); % k = 1 if i is in area j
end
for k = 1:narea(j)
  if temp(k) == i
    temp(k) = 0;
  end
end
temp = sort(temp);
temp = temp(2:narea(j));
area(j,1:narea(j)-1) = temp;
narea(j) = narea(j)-1;        
area = area(:,1:max(narea));  % remove trailing zero column
