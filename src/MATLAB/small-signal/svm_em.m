function [a] = svm_em(bus,line,opt)
% Syntax: [a] = svm_em(bus,line,opt)
%
% Purpose: generate linearized electromechanical model.
%          Note all loads are modeled as constant impedance. 
% 
% Input: bus - solved loadflow bus data
%        line - line data
%        opt - 0  to build full A matrix
%              1  to build only -inv(M)*K
%
% Output: a - system A matrix
%
% Files:
%
% See Also:

% Algorithm: 
%
% Calls:
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
% Date:     March 1992

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime 
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime 

jay = sqrt(-1);
if nargin == 2
  opt = 0;
end

[num_mach dummy] = size(mac_con);

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
Y_red       = red_ybus(bus,line);    % pre-fault
                                        % admittance matrix
% step 2: initialization
flag = 0;
f = mac_em(0,1,bus,flag);   
                   % all machine electro-mechanical model
mach_ref(1) = 0;

dmac_ang(1:num_mach,1) = zeros(num_mach,1);
dmac_spd(1:num_mach,1) = zeros(num_mach,1);

% step 3: perform perturbation 
if nargin <= 2 || opt == 0
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
 else
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
    a(1:num_mach,k) = dmac_spd(:,1)/pert;
    mac_ang(k,1) = nominal; 
  end
end  
