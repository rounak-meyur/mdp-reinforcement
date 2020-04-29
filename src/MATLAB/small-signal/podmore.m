function [n_bus,n_line,nmac_con]= ...
                      podmore(bus,line,area,nmach_a,basemva)
% Syntax   : [n_bus,n_line,nmac_con] = ...
%                     podmore(bus,line,area,nmach_a,basemva)
%
% Purpose  : To aggregate coherent machines using the 
%            Podmore method. A solved loadflow is required 
%            as input.
%
% Input    : bus      - bus data
%            line     - line data
%            area     - matrix of coherent machines
%            nmach_a  - vector of number of machines in each
%                         coherent area
%            basemva  - base mva (optional)
% Output   : n_bus  - new system bus data
%            n_line - new line data
%            nmac_con - aggregate generator data
%
% See also :
%
% Calls    :
%
% Call by  : 
%

% (c) Copyright 1991-4 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 1.0
% Author   : Pierre N. Accari and Joe H. Chow
% Date     : March 9, 1992
pst_var
jay=sqrt(-1);
if nargin == 4
  basemva = basmva;         % from global variable
end

t_bus  = bus;              % create new bus, line, 
t_line = line;             %   and machine data
t_mac  = mac_con;
[nrow,ncol] = size(mac_con);
if ncol <= 21
  t_mac = [t_mac ones(size(nrow,23-ncol))];
end

% set up internal bus index vector
nbus = length(t_bus(:,1));
bus_int = zeros(round(max(t_bus(:,1))),1);
for i = 1:nbus
  bus_int(t_bus(i,1)) = i;
end
% set up internal generator index vector
tot_mac  = length(mac_con(:,1)); % total number of machines
mac_int = zeros(round(max(mac_con(:,1))),1);
for i = 1:tot_mac
  mac_int(mac_con(i,1)) = i;
end

bus_vol  = t_bus(:,2);     % system bus voltages
bus_ang  = t_bus(:,3);     % system bus angles
num_area = length(area(:,1)); % get number of coherent areas

tot_cmac = 0.0;
bus_numb = t_bus(:,1);

nmac_con = [];
for r=1:num_area,             % cycle thru all areas
  if nmach_a(r) == 1          % single machine area
    nmac_con = [nmac_con; t_mac(mac_int(area(r,1)),:)];
    t_mac(mac_int(area(r,1)),1) = 0;  % set machine # to 0
   else                  % areas with more than 1 machine
    num_mach = nmach_a(r);   % # of coherent machines
    com_bus  = max(t_bus(:,1)) + 1;      % common bus 
       % number, make to one higher than largest bus number
    mac_list  = mac_int(area(r,1:nmach_a(r)))';   
                                 % coherent machine numbers
    new_mac  = area(r,1);
    bus_list1 = mac_con(mac_list,2);
    bus_list = bus_int(mac_con(mac_list,2));        
               % coherent machines bus numbers. It is okay 
               %   to have identical buses in bus_list.
    
    % search and eliminate lines connecting generator buses
    %   in the same coherent area
    nlines = length(t_line(:,1));
    line_list = [];
    for i = 1:nlines
      b1 = element(t_line(i,1),bus_list1);  % b1 and b2 = 1 
      b2 = element(t_line(i,2),bus_list1);  % => a match
      if b1+b2 == 1
        line_list = [line_list; i b1 b2];   % will use later
       elseif b1+b2 == 2  % line between two generator buses
        % replace line with injections
        bus1 = bus_int(t_line(i,1));
        v1 = bus(bus1,2)*exp(jay*bus(bus1,3)*pi/180);
        bus2 = bus_int(t_line(i,2));
        v2 = bus(bus2,2)*exp(jay*bus(bus2,3)*pi/180);
        [s1,s2] = line_pq(v1,v2,t_line(i,3),t_line(i,4), ...
                  t_line(i,5),t_line(i,6),t_line(i,7));
        t_bus(bus1,6) = t_bus(bus1,6) + real(s1);
        t_bus(bus1,7) = t_bus(bus1,7) + imag(s1);
        t_bus(bus2,6) = t_bus(bus2,6) + real(s2);
        t_bus(bus2,7) = t_bus(bus2,7) + imag(s2);
        t_line(i,1) = 0;  % 0 denotes line removed
      end
    end
    bus_type = 2;
    for ii = 1:num_mach
      if t_bus(bus_list(ii,1),10) == 1
        bus_type = 1;
      end
    end

    % inertia weighted aggregate machine internal voltage 
    %   and angle
    w_H = mac_con(mac_list,3).*mac_con(mac_list,16);
    % common bus (aggregate machine bus) voltage magnitude
    %   and angle
    mag_cbus = bus(bus_list,2)'*w_H/sum(w_H);   
    ang_cbus = bus(bus_list,3)'*w_H/sum(w_H);
    pfrac = mac_con(mac_list,22);
    qfrac = mac_con(mac_list,23);
    tot_pgen  = sum(t_bus(bus_list,4).*pfrac);
    tot_qgen  = sum(t_bus(bus_list,5).*qfrac);
    tot_pload = sum(t_bus(bus_list,6).*pfrac);
    tot_qload = sum(t_bus(bus_list,7).*qfrac);
    scale = (mag_cbus*ones(size(bus_list))./bus_vol(bus_list)).^2;
    tot_Gsh   = sum(t_bus(bus_list,8)./scale.*pfrac);
    tot_Bsh   = sum(t_bus(bus_list,9)./scale.*qfrac);

    %  add phase shifters to lines connecting to the 
    %   generator terminal buses
    nlines = length(line_list(:,1));
    for i = 1:nlines
      k = line_list(i,1);
      if line_list(i,2) == 1   % from bus matches
        if t_line(k,6) == 0
          t_line(k,6) = 1;
        end
        t_line(k,6) = t_line(k,6)*mag_cbus ...
                      /bus_vol(bus_int(t_line(k,1)));
        t_line(k,7) = t_line(k,7) + ang_cbus - ...
                      bus_ang(bus_int(t_line(k,1)));
        t_line(k,1) = com_bus;     % aggregate bus
       else                    % to bus matches
        temp = t_line(k,2);
        t_line(k,2) = t_line(k,1);
        if t_line(k,6) == 0
          t_line(k,6) = 1;
        end
        tap2 = t_line(k,6)^2;
        t_line(k,3) = t_line(k,3)*tap2;
        t_line(k,4) = t_line(k,4)*tap2;
        t_line(k,5) = t_line(k,5)/tap2;
        t_line(k,6) = (1/t_line(k,6))*mag_cbus ...
                      /bus_vol(bus_int(temp));
        t_line(k,7) = -t_line(k,7) + ang_cbus - ...
                      bus_ang(bus_int(temp));
        t_line(k,1) = com_bus;
      end
    end

    t_bus = [ t_bus; com_bus mag_cbus ang_cbus tot_pgen ...
     tot_qgen tot_pload tot_qload tot_Gsh tot_Bsh bus_type];
         
    % perform machine aggregation
    agg_mac = eqgen(t_mac,mac_list,100,com_bus,area(r,1));
    nmac_con = [nmac_con; agg_mac];
    t_mac(mac_list,1) = zeros(num_mach,1);  % set machine # 
                                            %  to 0
    t_bus(bus_list,1) = zeros(num_mach,1);  % set bus # to 0

  end
end

%  organize machine data
for i = 1:tot_mac
  if t_mac(i,1) ~= 0
    nmac_con = [nmac_con; t_mac(i,:)];
  end
end
%  organize bus data
num_bus = length(t_bus(:,1));
n_bus = [];
for i = 1:num_bus
  if t_bus(i,1) ~= 0
    n_bus = [n_bus; t_bus(i,:)];
  end
end
%  organize line data
num_line = length(t_line(:,1));
n_line = [];
for i = 1:num_line
  if t_line(i,1) ~= 0
    n_line = [n_line; t_line(i,:)];
  end
end
