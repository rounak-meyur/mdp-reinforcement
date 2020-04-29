% svm_16m.m
% m.file to generate the state variable models of the 16
%     machine, 68 bus system in LOTDYS report 
%     using the Matlab Power System Toolbox
jay = sqrt(-1);

% set up global variables
pst_var

% load input data from m.file
data16m
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 1 ;     % machine 1 is reference

% solve for loadflow
%   loadflow parameter
tol = 2e-13;   % tolerance for convergence
iter_max = 300; % maximum number of iterations
vmin = 0.5;   % voltage minimum
vmax = 1.5;   % voltage maximum
acc = 1.0;   % acceleration factor
[bus_sol,line_flw] = ...
      loadflow(bus,line,tol,iter_max,acc,'n',2);
bus = bus_sol;  % solved loadflow solution needed for initialization
save data16m bus line mac_con
load data16m
[num_mach, dummy] = size(mac_con);

%%
% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
Y_red       = red_ybus(bus,line);    % pre-fault admittance matrix

% step 2: initialization
flag = 0;
f = mac_em(0,1,bus,flag);   
                   % all machine electro-mechanical model
mach_ref(1) = 0;

dmac_ang(1:num_mach,1) = zeros(num_mach,1);
dmac_spd(1:num_mach,1) = zeros(num_mach,1);

%%
% step 3: perform perturbation 
for k = 1:num_mach
  dmac_ang(:,1) = zeros(num_mach,1);
  dmac_spd(:,1) = zeros(num_mach,1);
  pert = 0.001*mac_ang(k,1);   % one percent perturbation
  nominal = mac_ang(k,1);
  mac_ang(k,1) = mac_ang(k,1) + pert;
  % step 3a: network solution
  flag = 1;
  f = mac_em(0,1,bus,flag);
  psi = psi_re(:,1) + jay*psi_im(:,1);
  cur = Y_red*psi;
  cur_re(:,1) = real(cur); cur_im(:,1) = imag(cur);
  flag = 2;
  f = mac_em(0,1,bus,flag);
  a(1:num_mach,k) = dmac_ang(:,1)/pert;
  a(1+num_mach:2*num_mach,k) = dmac_spd(:,1)/pert;
  mac_ang(k,1) = nominal; 
end
for k = 1:num_mach
  dmac_ang(:,1) = zeros(num_mach,1);
  dmac_spd(:,1) = zeros(num_mach,1);
  pert = 0.001*mac_spd(k,1);   % one percent perturbation
  nominal = mac_spd(k,1);
  mac_spd(k,1) = mac_spd(k,1) + pert;
  % step 3a: network solution
  flag = 1;
  f = mac_em(0,1,bus,flag);
  psi = psi_re(:,1) + jay*psi_im(:,1);
  cur = Y_red*psi;
  cur_re(:,1) = real(cur); cur_im(:,1) = imag(cur);
  flag = 2;
  f = mac_em(0,1,bus,flag);
  a(1:num_mach,k+num_mach) = dmac_ang(:,1)/pert;
  a(1+num_mach:2*num_mach,k+num_mach) = dmac_spd(:,1)/pert;
  mac_spd(k,1) = nominal; 
end
MK = a(17:32,1:16);

save svm_16mk MK

