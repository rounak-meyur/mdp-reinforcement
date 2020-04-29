function [T] = sc_xform(area,nmach_a)
% Syntax   : [T] = sc_xform(area,nmach_a)
%
% Purpose  : To construct slow coherency transformation
%
% Input    : area     - matrix of coherent machines
%            nmach_a  - vector of number of machines in each
%                         coherent area
% Output   : T        - transformation matrix
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
% Version  : 1.0
% Author   : Joe H. Chow
% Date     : March 9, 1992

% set up internal generator index vector
tot_mac  = length(mac_con(:,1)); % total number of machines
mac_int = zeros(round(max(mac_con(:,1))),1);
for i = 1:tot_mac
  mac_int(mac_con(i,1)) = i;
end
num_area = length(area(:,1));

T1 = zeros(num_area,tot_mac); 
T2 = zeros(tot_mac-num_area,tot_mac);
T2_count = 0;
for r=1:num_area,             % cycle thru all areas
  if nmach_a(r) == 1          % single machine area
    T1(r,mac_int(area(r,1))) = 1;
   else                  % areas with more than 1 machine
    num_mach = nmach_a(r);   % # of coherent machines
    mac_list = mac_int(area(r,1:num_mach));
    w_H = mac_con(mac_list,3).*mac_con(mac_list,16);
    Ha = sum(w_H)
    for i = 1:num_mach
      T1(r,mac_list(i)) = mac_con(mac_list(i),3) ...
                   .*mac_con(mac_list(i),16)/Ha;
    end
    for i = 2:num_mach
      T2_count = T2_count+1;
      T2(T2_count,mac_list(1)) = -1;
      T2(T2_count,mac_list(i)) = 1;
    end
  end
end

T = [T1;T2];
