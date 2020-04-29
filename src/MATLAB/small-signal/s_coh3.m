function [n_bus,n_line,nmac_con] = ...
		       s_coh3(bus,line,area,nmach_a,basemva)
% Syntax   : [n_bus,n_line,nmac_con] = ...
%                      s_coh3(bus,line,area,nmach_a,basemva)
%
% Purpose  : To aggregate coherent machines using the
%            slow coherency method. In s_coh3, connection 
%            between the generator terminal buses are also
%            included for higher accuracy by a least squares
%            reconstruction of impedances. A solved loadflow
%            input is required.
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
% See also : podmore, i_agg 
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
% Author   : Joe H. Chow
% Date     : March 11, 1992
pst_var
jay=sqrt(-1);
if nargin == 4
  basemva = basmva;         % from global variable
end

n_bus  = bus;              % create new bus, line, 
n_line = line;             %   and machine data
n_mac  = mac_con;
[nrow,ncol] = size(mac_con);
if ncol <= 21
  n_mac = [n_mac ones(nrow,23-ncol)];
end

% set up internal bus index vector
nbus = length(n_bus(:,1));
bus_int = zeros(round(max(n_bus(:,1))),1);
for i = 1:nbus
  bus_int(n_bus(i,1)) = i;
end
% set up internal generator index vector
tot_mac  = length(mac_con(:,1)); % total number of machines
mac_int = zeros(round(max(mac_con(:,1))),1);
for i = 1:tot_mac
  mac_int(mac_con(i,1)) = i;
end

bus_vol  = n_bus(:,2);        % system bus voltages
bus_ang  = n_bus(:,3);        % system bus angles
num_area = length(area(:,1)); % number of coherent areas
nmac_con = [];
for r=1:num_area,             % cycle thru all areas
  if nmach_a(r) == 1          % single machine area
    nmac_con = [nmac_con; n_mac(mac_int(area(r,1)),:)];
    n_mac(mac_int(area(r,1)),1) = 0;  % set machine # to 0
   else                  % areas with more than 1 machine
    num_mach = nmach_a(r);   % # of coherent machines
    com_bus  = max(n_bus(:,1)) + 1;    % common bus 
       % number, make to one higher than largest bus number
    mac_list  = mac_int(area(r,1:nmach_a(r)))';   
				 % coherent machine numbers
    nlines   = length(n_line(:,1));    % # of lines
    bus_list1 = mac_con(mac_list,2);
    bus_list = bus_int(bus_list1); 
	       % coherent machines bus numbers. It is okay
	       %   to have identical buses in bus_list
    bus_type = 2;                     % generator bus
    for ii=1:num_mach,
      if n_bus(bus_list(ii,1),10) == 1,
	bus_type = 1;
      end
    end
    % inertia weighted aggregate machine internal voltage 
    %   and angle
    m = mac_con(mac_list,3).*mac_con(mac_list,16)/basemva;
    ma = sum(m);
    Pg = n_bus(bus_list,4).*mac_con(mac_list,22);
    Qg = n_bus(bus_list,5).*mac_con(mac_list,23);
    vb = n_bus(bus_list,2);
    angb = n_bus(bus_list,3)*pi/180;
    vbx = vb.*exp(jay*angb);
    xdp = basemva*mac_con(mac_list,7)./mac_con(mac_list,3);
    %  current and voltage computations include
    %  possibilities of multiple generator on the same bus
    int_cur = conj((Pg+jay*Qg)./vbx);
    int_vol = vbx + jay*xdp.*int_cur;
    % common bus (aggregate machine bus) voltage magnitude
    %   and angle
    mag_cbus = abs(int_vol)'*m/ma;   
    ang_cbus = angle(int_vol)'*m/ma;  
    vg = abs(int_vol);
    delta = angle(int_vol);
    % linearization
    K1 = -diag(vg.*vb.*cos(delta-angb)./xdp./m);
    K2 = [-diag(Pg./vb./m) -K1];
    K3 = -diag(int_vol./xdp);
    K4 = [ diag(exp(jay*angb)./(jay*xdp)) diag(vbx./xdp) ];
    %disp('stop at K1,K2,K3,K4')
    %disp('number of machines, machines')
    %r, mac_list
    %keyboard
    % Need to eliminate duplicate buses in bus_list due to multiple generators on the same bus
    % find the distinct buses in bus_list and store them in bus_listV
    bus_listV = []; bus_list2 = [];   % empty lists
    bus_listV(1,1) = bus_list(1,1);
    bus_list2(1,1) = bus_list1(1,1);
    for j = 2:num_mach
       add_status = element(bus_list(j,1),bus_listV);
       if add_status == 0
	  bus_listV = [bus_listV; bus_list(j,1)];
	  bus_list2 = [bus_list2; bus_list1(j,1)];
       end
    end
    num_dbus = length(bus_listV);       % number of distinct buses
    % set up contraction matrix
    TC = zeros(num_mach,num_dbus);
    for j = 1:num_mach
       for jj = 1:num_dbus
	  if bus_list(j,1) == bus_listV(jj,1)
	     TC(j,jj) = 1;
	  end
       end
    end

    % perform contraction
    K2 = K2 * [ TC 0*TC; 0*TC TC];
    K3 = TC' * K3;
    K4 = TC' * K4 * [ TC 0*TC; 0*TC TC];

    %keyboard

    % slow coherency transformation matrix
    T = [m'/ma; -ones(num_mach-1,1) eye(num_mach-1)]; 
			   % T = [C; G]

    % transform K1 and K2
    Tinv = inv(T);
    K1n = T*K1*Tinv;
    K2n = T*K2;
    K3n = K3*Tinv;
    %disp('see K1n,K2n,K3n')
    %keyboard
    % slow aggregate model
    invK1n22 = inv(K1n(2:num_mach,2:num_mach));
    K1a = ma*(K1n(1,1) - ...
	      K1n(1,2:num_mach)*invK1n22*K1n(2:num_mach,1));
    K2a = ma*(K2n(1,:) - ...
	      K1n(1,2:num_mach)*invK1n22*K2n(2:num_mach,:));
    K3a = K3n(:,1) - ...
	      K3n(:,2:num_mach)*invK1n22*K1n(2:num_mach,1);
    K4a = K4 - ...
	      K3n(:,2:num_mach)*invK1n22*K2n(2:num_mach,:);
    %disp('see K1a,K2a,K3a,K4a')    
    %keyboard
    % recover line parameters
    % aggregate generator
    beta1 = -K2a(1,1:num_dbus);
    beta2 = K2a(1,num_dbus+1:2*num_dbus);
    vb = n_bus(bus_listV,2);
    angb = n_bus(bus_listV,3)*pi/180;
    vbx = vb.*exp(jay*angb);
    phia = -atan2(beta1',beta2'./vb) + ...
		  ang_cbus*(1+0*vb) - angb;
    ra = sqrt( beta1'.^2 + (beta2'./vb).^2 );
    tapa = 1+0*vb;
    xdpa = mag_cbus./ra;
    v1 = mag_cbus*exp(jay*ang_cbus)*tapa;
    v2 = vb.*exp(jay*angb);
    [s1,s2] = line_pq(v1,v2,0*vb,xdpa,0*vb, ...
		      tapa,phia*180/pi);
    %disp('see xdpa, phia')
    %keyboard
    while(0) % 1 
    %  stop here to verify some calculations
    K3ax = -(v1.*exp(-jay*phia))./xdpa;
    %disp('K3a, K3ax, and (K3a-K3ax)')
    %[ K3a K3ax (K3a-K3ax)]
    disp('*****>'), pause
    disp('real power: Pg and real(s2) ')
    Pg,  real(s2)
    disp('*****>'), pause
    disp('reactive power: Qg and imag(s2) ')
    Qg, imag(s2)
    keyboard
    end  % while(0) of 1
    add_bus(1,1)   = com_bus;     % add common bus to
				      % bus data list
    add_bus(1,2)   = mag_cbus;    % common bus voltage
    add_bus(1,3)   = ang_cbus*180/pi; % angle
    add_bus(1,4:9) = zeros(1,6);
    add_bus(1,10)  = 3;           % bus type
    n_bus = [n_bus; add_bus];
    Pg = n_bus(bus_listV,4);
    Qg = n_bus(bus_listV,5);
    n_bus(bus_listV,4:5)= zeros(num_dbus,2); 
				      % remove Pgen & Qgen 
    for j = 1:num_dbus
      n_bus(bus_listV(j),6:7)= n_bus(bus_listV(j),6:7) - ...
	[Pg(j)+real(s2(j)) Qg(j)+imag(s2(j))]; % adjust load
    end
    n_bus(bus_list,10) = (1+0*bus_list)*3;   
				      % change load type
    % keyboard
    rl = nlines + num_dbus;
    n_line(nlines+1:rl,1) = (1+0*vb)*com_bus; 
					  % from bus
    n_line(nlines+1:rl,2) = n_bus(bus_listV,1); % to bus;
    n_line(nlines+1:rl,4) = xdpa;
    n_line(nlines+1:rl,6) = tapa; 
    n_line(nlines+1:rl,7) = phia*180/pi; 
    
    xdeq = 1/sum(ones(num_mach,1)./xdp);
    
    term_bus = max(n_bus(:,1))+1;     % new terminal bus # 
    com_vol = mag_cbus*exp(jay*ang_cbus);
    new_cur = conj(sum(s1)/com_vol);
    term_vol = com_vol -jay * xdeq * new_cur;

    add_bus(1,1)   = term_bus;
    add_bus(1,2)   = abs(term_vol);
    add_bus(1,3)   = angle(term_vol)*180/pi;
    add_bus(1,4)   = real(term_vol * conj(new_cur));
    add_bus(1,5)   = imag(term_vol * conj(new_cur)); 
    add_bus(1,6:9) = zeros(1,4); 
    add_bus(1,10)  = bus_type;
    n_bus = [n_bus; add_bus];

    n_line(rl+1,1)      = com_bus;
    n_line(rl+1,2)      = term_bus;
    n_line(rl+1,3)      = 0.0;
    n_line(rl+1,4)      = -xdeq;
    n_line(rl+1,5:7)    = zeros(1,3);

    %  recover lines between ex-generator buses
    K4_vb = K4a - [ diag(exp(jay*angb)./(jay*xdpa)) ...
		      diag(vbx./xdpa) ];
    % disp('see K4_vb?')
    % keyboard
    l_count = rl+1;
    for i = 1:num_dbus-1
      for j = i+1:num_dbus
	a = -K4_vb(i,j); b = -K4_vb(i,num_dbus+j);
	b_bar = b/(jay*vb(j));
	x = sqrt(a*b_bar);
	c = -K4_vb(j,i); d = -K4_vb(j,num_dbus+i);
	d_bar = d/(jay*vb(i));
	y = sqrt(c*d_bar);
	% assume phase shifter angle of 10 degrees
	% phi_ij = 10*pi/180;
	% assume phase shifter angle of 90 degrees
	phi_ij = pi/2;
	z1 = exp(jay*(angb(i)+angb(j))) ...
	     *(exp(jay*phi_ij)-exp(-jay*phi_ij)) ...
	     /(exp(jay*angb(i))*x-exp(jay*angb(j))*y);
	z2 = exp(jay*angb(j)) ...
	     /(x - exp(jay*(angb(j)+phi_ij))/z1);
	%disp('see a b b_bar x; c d d_bar y; z1 z2')
	%keyboard
	l_count = l_count + 1;
	n_line(l_count,:) = [bus_list2(i) bus_list2(j) ...
		    real(z1) imag(z1) 0 1. phi_ij*180/pi];
	l_count = l_count + 1;
	n_line(l_count,:) = [bus_list2(i) bus_list2(j) ...
		    real(z2) imag(z2) 0 1. 0];
	[si1,sj1]=line_pq(vbx(i),vbx(j), ...
		real(z1),imag(z1),0,1.,phi_ij*180/pi);
	[si2,sj2]=line_pq(vbx(i),vbx(j), ...
		real(z2),imag(z2),0,1.,0);
	%disp('see si1,sj1,si2,sj2')
        %keyboard
	%  match injections at buses i and j
	n_bus(bus_listV(i),6) = n_bus(bus_listV(i),6) ...
			       - real(si1+si2);
	n_bus(bus_listV(i),7) = n_bus(bus_listV(i),7) ...
			       - imag(si1+si2);
	n_bus(bus_listV(j),6) = n_bus(bus_listV(j),6) ...
			       - real(sj1+sj2);
	n_bus(bus_listV(j),7) = n_bus(bus_listV(j),7) ...
			       - imag(sj1+sj2);
        %disp('see bus injections')
        %keyboard
      end
      %disp('see n_bus, n_line')
      %keyboard
    end
% keyboard
    %  perform machine aggregation
    agg_mac = eqgen(n_mac,mac_list,basemva, ...
					term_bus,area(r,1));
    nmac_con = [nmac_con; agg_mac];
    n_mac(mac_list,1) = zeros(num_mach,1);  % set machine # 
					    % to zero
  end
end

%  organize machine data
for i = 1:tot_mac
  if n_mac(i,1) ~= 0
    nmac_con = [nmac_con; n_mac(i,:)];
  end
end

%disp('end of s_coh3.m')
%keyboard

