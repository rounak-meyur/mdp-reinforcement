function [area,nmach_a] = sign_coh(V_s)
% Syntax: [area,nmach_a] = sign_coh(V_s)
%
% Purpose: Find the coherent machine groups using the 
%          sign algorithm by Bo Eliasson.
%          (Ph.D. Thesis, Lund Institute, 1990)
%
% Input: V_s - eigenvectors corresponding to the slow
%              eigenvalues. The columns must be normalized
%              to one.
%
% Output: area - matrix of coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   coherent area
%
% See Also: V_slow, coh_map, L_group, ex_group, zabors, 
%           weaklink

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

[nmach,ns] = size(V_s);
%  stage 1 - construct sign eigenvector matrix 
for i = 1:nmach
  for j = 1:ns
    if V_s(i,j) > 0
       V_s(i,j) = 1;
      elseif V_s(i,j) < 0
       V_s(i,j) = -1;
    end
  end
end
%  compute grouping matrix
  G = V_s*V_s';
%  set up incidence matrix   
C = eye(ns);
for i = 1:nmach
  for j = 1:nmach
    if G(i,j) == ns
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
