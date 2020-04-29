function [n_mac]=eqgen(mac,mac_list,basemva,bus_num,mac_num)
% Syntax   : [n_mac]=eqgen(mac,mac_list,basemva, ...
%                          bus_num,mac_num)
%
% Purpose  : to construct an aggregate machine for 
%            a coherent area 
%
% Input    : mac      - machine data
%            mac_list - list of (internal) coherent machine 
%            basemva  - system base mva
%            bus_num  - aggregate machine bus number 
%            mac_num  - machine number to be assigned
%
% Output   : n_mac    - aggregate machine parameter
%
% See also : 
%
% Calls    : 
%
% Call by  :
%

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 2.0
% Author   : Joe H. Chow
% Date     : March 9, 1992
pst_var
num_mac = length(mac_list);         % # of coherent machines

%  new machine inertia on given mva base 
scales = mac(mac_list,16).*mac(mac_list,3)/basemva;
H_sum = sum(scales); 

n_mac(1,1) = mac_num;
n_mac(1,2) = bus_num;  n_mac(1,19) = bus_num;
n_mac(1,3) = basemva;
if min(mac(mac_list,4)) == 0;
    n_mac(1,4) = 0;
  else
    n_mac(1,4) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,4)/basemva);
end
if min(mac(mac_list,5)) == 0;
    n_mac(1,5) = 0;
  else
    n_mac(1,5) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,5)/basemva);
end
if min(mac(mac_list,6)) == 0;
    n_mac(1,6) = 0;
  else
    n_mac(1,6) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,6)/basemva);
end
if min(mac(mac_list,7)) == 0;
    n_mac(1,7) = 0;
  else
    n_mac(1,7) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,7)/basemva);
end
if min(mac(mac_list,8)) == 0;
    n_mac(1,8) = 0;
  else
    n_mac(1,8) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,8)/basemva);
end
n_mac(1,9) = scales'*mac(mac_list,9)/H_sum;
n_mac(1,10) = scales'*mac(mac_list,10)/H_sum;
if min(mac(mac_list,11)) == 0;
    n_mac(1,11) = 0;
  else
    n_mac(1,11) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,11)/basemva);
end
if min(mac(mac_list,12)) == 0;
    n_mac(1,12) = 0;
  else
    n_mac(1,12) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,12)/basemva);
end
if min(mac(mac_list,13)) == 0;
    n_mac(1,13) = 0;
  else
    n_mac(1,13) = 1/sum(mac(mac_list,3) ...
                 ./mac(mac_list,13)/basemva);
end
n_mac(1,14) = scales'*mac(mac_list,14)/H_sum;
n_mac(1,15) = scales'*mac(mac_list,15)/H_sum;
n_mac(1,16) = H_sum;
n_mac(1,17) = sum(mac(mac_list,17).*mac(mac_list,3)/basmva);
n_mac(1,18) = scales'*mac(mac_list,18)/H_sum;
n_mac(1,20) = scales'*mac(mac_list,20)/H_sum;
n_mac(1,21) = scales'*mac(mac_list,21)/H_sum;
n_mac(1,22) = 1;
n_mac(1,23) = 1;

