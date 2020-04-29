function [area,nmach_a,arear,nmach_ar] = ex_group(grouping)
% Syntax: [area,nmach_a,arear,nmach_ar] = ex_group(grouping)
%
% Purpose: Find the coherent machine groups using the slow
%          eigensubspace matrix. The number of coherent
%          groups is not restricted to be equal to the 
%          number of columns of V_s.
%          (Also called the tolerance-based slow-coherency
%           algorithm.)
%
% Input: grouping - output from coh_map, the difference of 
%                   the coherency map matrix and the
%                   coherency tolerance. 
%
% Output: area - matrix of tight coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   tight coherent area
%         arear - matrix of loose coherent groups
%         nmach_ar - vector of number of machines in each 
%                    loose coherent area
%
% See Also: V_slow, coh_map, L_group

% Algorithm: 
%
% Calls:
%
% Call By:

% (c) Copyright 1991-4 Joe H. Chow - All Rights Reserved

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

[nmach,dum] = size(grouping);
%  stage 1 - loose coherent groups
%  find the coherent machines for each machine based on 
%  the minimum tolerance
n_coh = ones(nmach,1);
coh_area = [1:1:nmach]';
coh_label = coh_area;
for i = 1:nmach
  for j = i+1:nmach
    if grouping(i,j) > 0 % machines i & j are coherent
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
arear = coh_area(select,:);
nmach_ar = n_coh(select);
% sort machines in ascending order
for i = 1:nr
  temp = arear(i,1:nmach_ar(i));
  temp = sort(temp);
  arear(i,1:nmach_ar(i)) = temp;
end
% keyboard
%
% stage 2 - tight coherent groups
% fine tune the coherent groups
n = 0;
disp('number of loose-coherent areas')
nr 
for i = 1:nr
%  disp('loose area')
%  i
  if nmach_ar(i) == 1 % single machine area, no splitting
    n = n+1;
    area(n,1) = arear(i,1);
    nmach_a(n) = nmach_ar(i);
   else
    temp = arear(i,1:nmach_ar(i));
    % extract appropriate part of coh_map
    map = grouping(temp',:);
    map = map(:,temp');
    % check if the coherent group is tight
    row_sum = sum(map - diag(diag(map)));
    if min(row_sum) > 0
      n = n+1;
      area(n,1:nmach_ar(i)) = temp;
      nmach_a(n) = nmach_ar(i);
     else
      iarea = 1;      % number of machine areas
      nmach = nmach_ar(i);
      narea_tem = nmach;  % number of machines in each area
      area_tem = [1:1:nmach]; % machines in each area
         % machine numbers are reassigned from 1
      cost_old = 0;
      term = 0;
      while term == 0
        % set off-diagonal blocks to zero
        if iarea == 1 
          map_diag = map;
         else
          map_diag = diagb(map,area_tem,narea_tem);
        end
        % find smallest row sum in map_diag, ie, identify
        %   least coherent machine
        y = sum(map_diag);
        for k = 1:iarea
           jk = area_tem(k,1:narea_tem(k));
           y(jk) = y(jk)/narea_tem(k);
        end
        [dum,next_mac] = min(y);
        % take next_mac from its area
        [areay,nareay] = ...
                    rem_area(area_tem,narea_tem,next_mac);
        cost = 0;
%        disp('iarea')
%        iarea
        for k = 1:iarea+1
          [areax,nareax]=add_area(areay,nareay,next_mac,k);
          % zero out diagonal blocks
          map_off =  offdb(map,areax,nareax);
          new_cost = sum(sum(map_off));
          if new_cost < cost
            cost = new_cost;
            area_tem = areax;
            narea_tem = nareax;
          end
        end
        iarea = length(narea_tem);
        if cost < cost_old
          cost_old = cost;
         else
          term = 1;     % terminate while loop 
          n_add = size(narea_tem);
          for k = 1:n_add
            n = n+1;  
            area(n,1:narea_tem(k)) = ...
                       temp(area_tem(k,1:narea_tem(k)));
            nmach_a(n) = narea_tem(k);
          end
        end
      end
    end
  end

end

disp('number of tight-coherent areas')
length(nmach_a) 
