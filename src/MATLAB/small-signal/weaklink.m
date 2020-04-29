function [area,nmach_a,areal,nmach_al,list,S] = ...
                                          weaklink(MK,tol)
% Syntax: [area,nmach_a,areal,nmach_al,list,S] = ...
%                                         weaklink(MK,tol)
%
% Purpose: Find the coherent machine groups using the 
%          weak link method by Nath, Lamba & Rao.
%          (IEEE Trans. PAS, vol. 104, pp. 1443-9, 1985.)
%
% Input: MK - the linearized -inv(M)*K matrix, where M is 
%             the machine inertia, and K is the 
%             synchronizing coefficients.
%        tol - tolerance for elimination of weak coupling
%
% Output: area - matrix of strong coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   strong coherent area
%         areal - matrix of weak coherent groups of machines
%         nmach_al - vector of number of machines in each 
%                    weak coherent area
%         list - machine order based on coupling factor
%         S    - coupling factors
%
% See Also: V_slow, coh_map, L_group, ex_group, zabors

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

%  stage 1 - sort the machines according to S
[nmach,dum] = size(MK);
%  normalize MK with respect to its largest in magnitude entry
MKabs = abs(MK);
STOT = sum(sum(MKabs));
%  select 1st machine
x = diag(MKabs);
[x,ind] = sort(x);
list1 = ind(1);
list2 = ind(2:nmach);
SNUM = sum(MKabs(:,ind(1))) + sum(MKabs(ind(1),:)) ...
       - 2*MKabs(ind(1),ind(1));
SDEN = STOT - SNUM;
S = SNUM/SDEN;

%  sort thru other machines
for i = 2:nmach-1
  nj = length(list2);
  Smin = STOT;
  for j = 1:nj
    tlist1 = [list1; list2(j)];
    tlist2 = list2(1:nj-1);
    tlist2(j) = list2(nj);
    SNUM = 0;
    for k = 1:length(tlist1)
      SNUM = SNUM + sum(MKabs(tlist1(k),tlist2)) + ...
                    sum(MKabs(tlist2,tlist1(k)));
    end
    if SNUM < Smin
      Smin = SNUM;
      indj = j;
    end
  end
  list1 = [list1; list2(indj)];
  tmp = list2(nj);
  list2 = list2(1:nj-1);
  if indj ~= nj
    list2(indj) = tmp;
  end
  SDEN = STOT - Smin;
  S = [S; SNUM/SDEN];
end
list1 = [list1; list2];
S = [0; S; 0];
plot(S), title('Coupling Factors')
keyboard

% find strong cokherent areas
for i = 1:nmach
  dS(i) = S(i+1) - S(i);
end

for i = 1:nmach-1
  d2S(i) = dS(i+1) - dS(i);
end

area = list1(1);
nmach_a = 1;
num_area = 1;
for i = 1:nmach-1
  if d2S(i) > 0   % generate new area
    num_area = num_area + 1;
    area(num_area,1) = list1(i+1);
    nmach_a = [nmach_a; 1];
   else           % add machine to last area
    nmach_a(num_area) = nmach_a(num_area) + 1;
    area(num_area,nmach_a(num_area)) = list1(i+1);
  end
end

if nargin == 1
   areal = list1;
   nmach_al = S;
  else
   % find weak coherent areas
   areal = list1(1);
   nmach_al = 1;
   num_area = 1;
   for i = 1:nmach-1
     if d2S(i) > tol   % generate new area
       num_area = num_area + 1;
       areal(num_area,1) = list1(i+1);
       nmach_al = [nmach_al; 1];
      else           % add machine to last area
       nmach_al(num_area) = nmach_al(num_area) + 1;
       areal(num_area,nmach_al(num_area)) = list1(i+1);
     end
   end
   bar(d2S), title('Grouping Bar chart')
   keyboard
   list = list1;
end
