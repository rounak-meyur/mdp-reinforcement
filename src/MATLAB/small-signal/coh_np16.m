% coh_np16.m
% m.file to test the coherency grouping methods for the 
%     NPCC 16 machine, 68 bus system in LOTDYS report 
%     using the Matlab Power System Toolbox
clear all

jay = sqrt(-1);

% set up global variables
pst_var

disp('Coherency demo program')
disp('The system is the NPCC 16 machine, 68-bus system')
disp('Solving loadflow')

% load input data from m.file
data16m
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 1 ;     % machine 1 is reference

% solve for loadflow
%   loadflow parameter
tol =1e-12;   % tolerance for convergence
iter_max = 30; % maximum number of iterations
vmin = 0.5;   % voltage minimum
vmax = 1.5;   % voltage maximum
acc = 1.0;   % acceleration factor
[bus_sol,line_flw] = ...
       loadflow(bus,line,tol,iter_max,acc,'n',2);
bus = bus_sol;  % solved loadflow solution needed for
                 % initialization
% save data16m bus line mac_con
% load data16m

q = 0;
while( q == 0)
  disp('enter 1 for slow coherency method')
  disp('      2 for tight coherency method')
  disp('      3 for Zaborszky clustering method')
  disp('      4 for Rao weak link method')
  disp('      5 for Eliasson sign coherency method')
  disp('      0 to quit')
  sel = input('enter selection >> ');

  if sel == 1
     %  slow coherency grouping
     disp('5 slowest eigenvalues will be used')
     [area,nmach_a] = coherent(bus,line,1,5,0);
     disp('coherent areas')
     area
   elseif sel == 2 
     %  tight slow coherency grouping
     disp('5 slowest eigenvalues and tolerance = 0.95 ')
     disp('  will be used')
     [area,nmach_a,areal,nmach_al]= ... 
                             coherent(bus,line,2,5,.95);
     disp('coherent areas')
     area
   elseif sel == 3
     %  Zaborszky's clustering technique
     disp('alpha = 0.5 will be used')
     [area,nmach_a] = coherent(bus,line,3,0,0.5);
     disp('coherent areas')
     area
   elseif sel == 4
     %  weak link technique
     disp('tolerance = 0.05')
     [area,nmach_a,areal,nmach_al]= ...
                             coherent(bus,line,4,0,.05);
     disp('coherent areas')
     area
   elseif sel == 5
     %  Eliasson technique
     disp('5 slowest eigenvalues will be used')
     [area,nmach_a]=coherent(bus,line,5,5,0);
     disp('coherent areas')
     area
   elseif sel == 0
     q = 1;
   else
     disp('invalid selection')
  end
end

