% agg_np48.m   Joe Chow 7/92
% m file to generate a reduced order power system model 
%   for the NPCC 48 machine system for the Medway 
%   disturbance.
clear
jay = sqrt(-1);

pst_var                      % set up global variables

disp('Aggregation demo program')
disp('The system is the NPCC 48 machine system. ')
disp('Reading loadflow')

datanp48                     % load input data from m.file
basrad  = 2*pi*60;           % system frequency is 60 Hz
basmva  = 100;               % 100 MVA base
load datanp48                % solved loadflow file

disp(' ')
disp('Reduced order model is for Medway disturbance.')
disp('  New England area kept in full detail.')
disp('  Reduced model will have 24 machines.')
disp('  If all the areas are to be aggregated, just comment')
disp('  out two lines of code in agg_np48.m .')

% the following coherent areas are the 17 area partition 
%   obtained from ex_group

area = [ ...
  3  4  5  6  7  8  0  0  0;
  1  2  9  0  0  0  0  0  0;
 10  0  0  0  0  0  0  0  0;
 11 12  0  0  0  0  0  0  0;
 13 14 24 25 26  0  0  0  0;
 15 16 17 18 19 20 21 22 23;
 27 28 29 30  0  0  0  0  0;
 31  0  0  0  0  0  0  0  0;
 32 37 38 40 42  0  0  0  0;
 33  0  0  0  0  0  0  0  0;
 34 35  0  0  0  0  0  0  0;
 36  0  0  0  0  0  0  0  0;
 39  0  0  0  0  0  0  0  0;
 41  0  0  0  0  0  0  0  0;
 43 44 45 46  0  0  0  0  0;
 47  0  0  0  0  0  0  0  0;
 48  0  0  0  0  0  0  0  0 ];
narea = [ 6 3 1 2 5 9 4 1 5 1 2 1 1 1 4 1 1 ];

% take New England out of area; comment the following two lines 
%   out if all areas are to be aggregated.
area = [area(3:17,:)];
narea = [narea(3:17)];
disp(' ')
disp('enter 1 for Podmore aggregation method')
disp('      2 for inertial aggregation method')
disp('      3 for slow coherency method')
sel = input('enter selection >> ');

if sel == 1 
    % perform podmore aggregation
    [rbus,rline,rmac]=podmore(bus,line,area,narea,100);
  elseif sel == 2
    % perform inertia aggregation
    [rbus,rline,rmac]=i_agg(bus,line,area,narea,100);
  elseif sel == 3
    % perform slow coherency aggregation
    [rbus,rline,rmac]=s_coh3(bus,line,area,narea,100);
  else 
    error('invalid selection')
end
disp(' ')
disp('Reduced order model constructed')
disp('  checking loadflow for reduced order model')
% check loadflow 
  tol      = 1e-9;           % tolerance for convergence
  iter_max = 30;             % maximum number of iterations
  vmin     = 0.5;            % voltage minimum
  vmax     = 1.5;            % voltage maximum
  acc      = 1.0;            % acceleration factor
  [bus_sol,line_flw] = loadflow(rbus,rline,tol, ...
                             iter_max,vmin,vmax,acc,'n',2); 
rbus = bus_sol;

bus = bus_sol;
line = rline;
mac_con = rmac;

disp(' ')
disp('Generating linearized electromechanical model')
a = svm_em(bus,line,0);
disp(' ')
disp('Computing eigenvalues, stored in the array modes')
modes = sort(eig(a))

disp(' ')
disp('Eliminating the odd number buses in New England')

bus_list = [1:2:29]';
flag = 0; % eliminate buses in bus_list
[bus_red,line_red]=reduce(bus,line,bus_list,flag);

disp(' ')
disp('Checking loadflow for reduced network')
[bus_sol,line_flow] = loadflow(bus_red,line_red, ...
                                   1e-10,5,.5,1.5,1,'n',2);

