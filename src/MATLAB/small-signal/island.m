function [bus_isl,line_isl,bbus,mac_isl,exc_isl] = ...
                                  island(bus,line,bus_list)
% Syntax: [bus_isl,line_isl,bbus,mac_isl,exc_isl]
%                                =island(bus,line,bus_list)
%
% Purpose: Extract a subsystem from a power system
%
% Input: bus      - solved loadflow bus data
%        line     - line data
%        bus_list - buses in the subsystem
%
% Output: bus_isl - bus data for buses in the subsystem
%         line_isl - line data for subsystem
%         bbus - boundary buses
%         mac_isl - machine data for the subsystem
%         exc_isl - exciter data for the subsystem
%
% See Also: 

% Algorithm: 
%
% Calls: element
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
% Date:     August 1991
pst_var
% get bus data
[nbus dum] = size(bus);
bus_isl = [];
for i = 1:nbus
  if element(bus(i,1),bus_list) == 1
    bus_isl = [bus_isl ; bus(i,:)];
  end
end

% get line data and boundary bus data
[nline dum] = size(line);
line_isl = [];
bbus = [];
for i = 1:nline
  if element(line(i,1),bus_list) == 1 & ...
     element(line(i,2),bus_list) == 1
      line_isl = [line_isl ; line(i,:)];
    elseif element(line(i,1),bus_list) == 1
      if isempty(bbus) == 1
         bbus = line(i,1);
        else
         if element(line(i,1),bbus) == 0
           bbus = [bbus; line(i,1)];
         end
      end
    elseif element(line(i,2),bus_list) == 1
      if isempty(bbus) == 1
         bbus = line(i,2);
        else
         if element(line(i,2),bbus) == 0
           bbus = [bbus; line(i,2)];
         end
      end
  end
end
bbus = sort(bbus);

% get machine data 
mac_isl = [];
mac_list = [];
if exist('mac_con') == 1
  if isempty(mac_con) == 0
    [nmach dum] = size(mac_con);
    for i = 1:nmach
      if element(round(mac_con(i,2)),bus_list) == 1
        mac_isl = [mac_isl ; mac_con(i,:)];
        mac_list = [mac_list ; mac_con(i,1)];
      end
    end
  end
end

% get exciter data 
exc_isl = [];
if exist('exc_con') == 1
  if isempty(exc_con) == 0
    [nexc dum] = size(exc_con);
    for i = 1:nexc
      if element(round(exc_con(i,2)),mac_list) == 1
        exc_isl = [exc_isl ; exc_con(i,:)];
      end
    end
  end
end
