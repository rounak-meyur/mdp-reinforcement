%  np48map.m
%    m file to generate coherency map

clear 
pst_var

% load input data from m.file
datanp48
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 1 ;     % machine 1 is reference

%while(0)
% solve for loadflow
%   loadflow parameter
tol = 2e-13;   % tolerance for convergence
iter_max = 30; % maximum number of iterations
vmin = 0.5;   % voltage minimum
vmax = 1.5;   % voltage maximum
acc = 1.0;   % acceleration factor
[bus_sol,line_flw] = ...
       loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
bus = bus_sol;  % solved loadflow solution needed for
                 % initialization
save datanp48 bus line mac_con
%end % while(0)

load datanp48

% generate -inv(M)*K matrix
MK = svm_em(bus,line,1);
% extract first ns eigenvalues
ns = 9;
[eig_s,V_s] = V_slow(MK,ns);

% [area,nmach_a,L_x] = L_group(V_s);
tol = 0.95
[nrow ncol] = size(V_s);
[grouping,c_map] = coh_map(V_s,tol);
save c_map c_map
m = [30 60];  % orientation for 3D plot
mesh(20*c_map)
xlabel('Machine number')
ylabel('Machine number')
zlabel('Coherency measure')
%ylabel('Values of C_m')
%gtext('Columns of C_m')
%gtext('Rows of C_m')
%gtext('(1,1)')

disp(' ')
disp('enter any key to continue'), pause

%  ask for 3D mode shape display

mode = 1;
while(mode)
  disp('3D mode shape display')
  mode = input('Enter mode number to see or 0 to stop > ')
  if mode ~= 0
    mode_map = bus_map;
    for i = 1:nrow
      j = mac_con(i,2);
      for k = 1:length(bus_xy(:,1))
        if bus_xy(k,1) == j 
          mode_map(bus_xy(k,2),bus_xy(k,3)) = ... 
                sign(mode)*V_s(i,abs(mode));
% new code
          mode_map(bus_xy(k,2)-1,bus_xy(k,3)-1) = ... 
                sign(mode)*V_s(i,abs(mode));
          mode_map(bus_xy(k,2),bus_xy(k,3)-1) = ... 
                sign(mode)*V_s(i,abs(mode));
          mode_map(bus_xy(k,2)-1,bus_xy(k,3)) = ... 
                sign(mode)*V_s(i,abs(mode));
% end new code
        end
      end
    end
save mode_map mode_map
    disp('type < mesh(mode_map) > to plot')
    keyboard
    xlabel('x grid')
    ylabel('y grid')
    zlabel('Mode shape')
  end
end

