function [area,nmach_a] = zabors(MK,tol)
% Syntax: [area,nmach_a] = zabors(MK,tol)
%
% Purpose: Find the coherent machine groups using the 
%          clustering algorithm by John Zaborszky.
%          (IEEE Trans. CAS, vol. 29, pp. 747-757, 1982.)
%
% Input: MK - the linearized -inv(M)*K matrix, where M is 
%             the machine inertia, and K is the 
%             synchronizing coefficients.
%        tol - tolerance for elimination of weak coupling
%
% Output: area - matrix of coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   coherent area
%
% See Also: V_slow, coh_map, L_group, ex_group

% Algorithm: 
%
% Calls:
%
% Call By:

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1992

%  stage 1 - construct reduced incidence matrix
[nmach,dum] = size(MK);
%  normalize MK with respect to its largest in magnitude entry
MKabs = abs(MK);
MKmax = max(max(MKabs));
MKabs = MKabs/MKmax;
%  go thru each column and eliminate small elements
for i = 1:nmach
  rowi = MKabs(i,:);
  [rowi,ind] = sort(rowi);
  alpha = tol;
  for j = 1:nmach
    alpha = alpha - rowi(j);
    if alpha < 0
      break
    end
    MKabs(i,ind(j)) = 0;
  end
end
%  set up reduced incidence matrix   
C = eye(size(MK));
for i = 1:nmach
  for j = 1:nmach
    if MKabs(i,j) > 0
      C(i,j) = 1;
      C(j,i) = 1;
    end
  end
end

%  stage 2 - find coherent groups
n_coh = ones(nmach,1);
coh_area = [1:1:nmach]';
coh_label = coh_area;
for i = 1:nmach
  for j = i+1:nmach
    if C(i,j) > 0 % machines i & j are coherent
      if coh_label(i) ~= coh_label(j) 
        k1 = coh_label(i);
        k2 = coh_label(j);
        k_small = min(k1,k2);
        k_large = max(k1,k2);
        coh_label(coh_area(k_large,1:n_coh(k_large))') = ...
                            k_small*ones(n_coh(k_large),1);
        temp = [ coh_area(k_small,1:n_coh(k_small)) ...
                 coh_area(k_large,1:n_coh(k_large)) ];
        n_coh(k_small)= n_coh(k_small)+n_coh(k_large);
        coh_area(k_small,1:n_coh(k_small)) = temp;
        coh_area(k_large,1:n_coh(k_large)) = ...
                             zeros(1,n_coh(k_large));
        n_coh(k_large) = 0;
      end
    end
  end
end
select = [];
for i = 1:nmach
  if n_coh(i) ~= 0
    select = [select ; i];
  end
end
nr = length(select);
area = coh_area(select,:);
nmach_a = n_coh(select);
% sort machines in ascending order
for i = 1:nr
  temp = area(i,1:nmach_a(i));
  temp = sort(temp);
  area(i,1:nmach_a(i)) = temp;
end
