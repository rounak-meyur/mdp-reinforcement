function [bus_red,line_red] = reduce(bus,line, ...
                                          bus_list,flag,tol)
% Syntax : [bus_red,line_red] = ...
%                         reduce(bus,line,bus_list,flag,tol)
% Purpose: performs a network reduction on the system. A 
%          list of buses contained in bus_list will either 
%          be retained or eliminated.
%
% Input  : bus      - system bus data
%          line     - system line data
%          bus_list - list of buses to be retained or eliminated
%          flag     - 1 to retain buses in bus_list
%                     0 to eliminate buses in bus_list
%          tol      - impedance tolerance
%
% Output : bus_red  - reduced bus list
%          line_red - reduced line list
%
% See also :
%
% Calls    :
%
% Call by  :
%

% (c) Copyright 1992 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 1.0
% Author   : Pierre N. Accari and Joe H. Chow
% Date     : 4 March, 1992
pst_var
jay = sqrt(-1);

%len_bus  = length(bus);  % this statement is wrong
len_bus  = length(bus(:,1));  
[m n]    = size(bus_list);
bus_temp = zeros(len_bus,2);

if nargin == 3,        % if there are three input arguments
   flag = 1;           % set flag to retain
   tol  = 1e-8;        % set tolerance to value

elseif nargin == 4,    % if there are four input arguments
   tol  = 1e-8;        % set tolerance to value
end

for i=1:len_bus,       % set up internal bus numbering
   bus_int(round(bus(i,1))) = i;
end

bus_temp = [1:1:len_bus]';
if flag == 1,   % buses to be retained. 
   bus_temp(:,2) = ones(len_bus,1); % 1 is eliminate
   bus_temp(bus_int(bus_list),2) = zeros(size(bus_list));
  elseif flag == 0, % buses to be eliminated
   bus_temp(:,2) = zeros(len_bus,1); % 1 is eliminate
   bus_temp(bus_int(bus_list),2) = ones(size(bus_list));
end

   cnt_temp = 0;
   for i=1:len_bus,   % retain all generator buses
       if ( bus(i,10) == 1 | bus(i,10) == 2 ),
          cnt_temp = cnt_temp + 1;
          bus_temp(i,2) = 0;
       end
   end

   [x I] = sort(bus_temp(:,2));   % sort bus list so that
     % all buses to be eliminated are at the end of the list
   bus_temp = [bus_temp(I,1) bus_temp(I,2)];

   s_blist = 0;             % number of retained buses
   l_blist = length(bus_temp);
   for i=1:l_blist
      if bus_temp(i,2) == 0,
         s_blist = s_blist +1;
      end
   end
                                                            % eliminated buses
   bus_list = bus_temp(1:s_blist,1);

bus_ret = bus(bus_list,:);
bus_eli = bus(bus_temp(s_blist+1:len_bus,1),:);
bus_ord = [ bus_ret ; bus_eli];


max_line = length(line(:,1));

l_blist = length(bus_list);
ret_cnt  = 1;
eli_cnt  = 1;

for i=1:max_line,   % check line data. If both from_bus and 
   ff_bus   = 0;    % to_bus arein the list of buses to be 
   tf_bus   = 0;    % retained, set flag to retain so that 
   from_bus = line(i,1); % corresponding line will be 
   to_bus   = line(i,2); % retained for the study.
   for j=1:l_blist,
      if bus_int(from_bus) == bus_list(j),
         ff_bus = 1;
      end
      if bus_int(to_bus) == bus_list(j),
         tf_bus = 1;
      end
   end

   if (ff_bus == 1) & (tf_bus == 1),    % set up retained 
      line_ret(ret_cnt,:) = line(i,:);  % line data
      ret_cnt = ret_cnt + 1;
     else                      % set up eliminated line data
      line_eli(eli_cnt,:) = line(i,:);
      eli_cnt = eli_cnt + 1;
   end
end

line_ord = [ line_ret ; line_eli ];

ij = 0;                         % line count

[Y] = Ybus(bus_ord,line_eli);   %  create admittance matrix
[mb nb] = size(bus_list);
[my ny] = size(Y);
%  add loads of buses to eliminated in Y
for i = mb+1:ny
  Y(i,i) = Y(i,i) + (bus_ord(i,6)-bus_ord(i,4) ...
           -jay*(bus_ord(i,7)-bus_ord(i,5)))/bus_ord(i,2)^2;
end
Y_11 = Y(1:mb,1:mb);      % patition matrix into submatrices
Y_12 = Y(1:mb,mb+1:ny);
Y_21 = Y(mb+1:my,1:mb);
Y_22 = Y(mb+1:my,mb+1:ny);
Y_red = Y_11 - Y_12 * inv(Y_22) * Y_21; % reduce network
bus_red = bus_ret;

ij = 0;
for i=1:mb        % reconstruct system from reduced matrix
  for j=i+1:mb 
    if abs(Y_red(i,j)) >= tol
      if abs(abs(Y_red(i,j)) - abs(Y_red(j,i))) <= tol,
%     if abs(Y_red(i,j)) == abs(Y_red(j,i)),
      % only one line needed
        if abs(real(Y_red(i,j))-real(Y_red(j,i)))<= tol & ...
           abs(imag(Y_red(i,j))-imag(Y_red(j,i))) <= tol 
           % no phase shifter needed
            ij = ij + 1;
            y = - ( Y_red(i,j) + Y_red(j,i) ) / 2;
            z = 1 / y;
            line_rec(ij,1) = bus_red(i,1);
            line_rec(ij,2) = bus_red(j,1);
            line_rec(ij,3) = real(z);
            line_rec(ij,4) = imag(z);
            line_rec(ij,5) = 0;
            line_rec(ij,6) = 0;
            line_rec(ij,7) = 0;
            Y_red(i,i) = Y_red(i,i) - y;
            Y_red(j,j) = Y_red(j,j) - y;

         else

            ij = ij + 1;
            alpha = Y_red(i,j)/Y_red(j,i);
            theta = atan2(imag(alpha),real(alpha))/2;
            y = -Y_red(j,i)*exp(jay*theta);
            z = 1 / y;
            line_rec(ij,1) = bus_red(i,1);
            line_rec(ij,2) = bus_red(j,1);
            line_rec(ij,3) = real(z);
            line_rec(ij,4) = imag(z);
            line_rec(ij,5) = 0;
            line_rec(ij,6) = 1;
            line_rec(ij,7) = theta*180/pi;
            Y_red(i,i) = Y_red(i,i) - y;
            Y_red(j,j) = Y_red(j,j) - y;
          
        end  

       else
       % two line reconstruction
            ij = ij + 1;
            y_1 = - ( Y_red(i,j) + Y_red(j,i) ) / 2;
            z_1 = 1 / y_1;
            line_rec(ij,1) = bus_red(i,1);
            line_rec(ij,2) = bus_red(j,1);
            line_rec(ij,3) = real(z_1);
            line_rec(ij,4) = imag(z_1);
            line_rec(ij,5) = 0;
            line_rec(ij,6) = 0;
            line_rec(ij,7) = 0;
            Y_bar_1 = Y_red(i,j) + y_1;
            Y_bar_2 = Y_red(j,i) + y_1;
            alpha = Y_bar_1/Y_bar_2;
            phi = atan2(imag(alpha),real(alpha))/2;
            y_2     = - Y_bar_1 * exp(-jay * phi);
            z_2     = 1 / y_2;
            ij = ij + 1;
            line_rec(ij,1) = bus_red(i,1);
            line_rec(ij,2) = bus_red(j,1);
            line_rec(ij,3) = real(z_2);
            line_rec(ij,4) = imag(z_2);
            line_rec(ij,5) = 0;
            line_rec(ij,6) = 1;
            line_rec(ij,7) = phi*180/pi;
            Y_bar_11 = y_1 + y_2;
            Y_red(i,i) = Y_red(i,i) - Y_bar_11;
            Y_red(j,j) = Y_red(j,j) - Y_bar_11;

      end

    end
  end
  %  recover bus shunts
  bus_red(i,8) = real(Y_red(i,i));
  bus_red(i,9) = imag(Y_red(i,i));
end
line_red = [line_ret; line_rec];

