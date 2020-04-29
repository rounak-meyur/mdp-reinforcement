function [area,nmach_a,L] = L_group(V_s)
% Syntax: [area,nmach_a,L] = L_group(V_s)
%
% Purpose: Find the coherent machine groups using the slow
%          eigensubspace matrix. The number of coherent
%          groups is equal to the number of columns of V_s
%
% Input: V_s - eigenvectors corresponding to the slow
%              eigenvalues. The columns must be normalized
%              to one.
%
% Output: area - matrix of coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   area
%         L - matrix L of non-reference machines
%
% See Also: V_slow, coh_map

% Algorithm: See J. H. Chow, "Time-scale modeling of dynamic
%            networks with applications to power systems,"
%            Chapter 5.
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

[nmach,narea] = size(V_s);
V = [zeros(1,narea+1); zeros(nmach,1) V_s];
mach_ord = [1:1:nmach]';
for i = 1:narea   % find reference machine for each area
  V = V(2:nmach-i+2,2:narea-i+2);
  [V_max,row_ind] = max(abs(V)); % find max in each column
                                 % for reference machine
  [V_max,col_ind] = max(V_max);  % find overall max
  % rearrange the machine ordering
  new_ref = row_ind(col_ind);
  temp = mach_ord(i);
  mach_ord(i) = mach_ord(new_ref+i-1);
  mach_ord(new_ref+i-1) = temp;
  % permute the V matrix
  temp = V(1,:);          % 1st do column permutation
  V(1,:) = V(new_ref,:); 
  V(new_ref,:) = temp;
  temp = V(:,1);          % next do row permutation
  V(:,1) = V(:,col_ind);
  V(:,col_ind) = temp;
  % perform Gaussian elimination
  v = V(:,1)/V(1,1);
  X = v*V(1,:);
  V = V - X;
end
% permute the V_s matrix according to mach_ord
temp1 = mach_ord(1:narea);
temp2 = mach_ord(narea+1:nmach);
temp2 = sort(temp2);
mach_ord = [temp1; temp2]; % non-reference machines in
                           % ascending order
V_new = V_s(mach_ord,:);
% partition the V_new matrix
V1 = V_new(1:narea,:); V2 = V_new(narea+1:nmach,:);
L = V2/V1; % V2*inv(V1);
[dum,k] = max(L'); % k is area assignment vector
nmach_a = ones(1,narea);
area(:,1) = mach_ord(1:narea);
%  assign non-reference machines to areas
for i = 1:nmach-narea
  nmach_a(k(i)) = nmach_a(k(i)) + 1;
  area(k(i),nmach_a(k(i))) = mach_ord(narea+i);
end

