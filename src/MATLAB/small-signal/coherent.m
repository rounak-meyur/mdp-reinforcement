function [area,nmach_a,areal,nmach_al] = ...
                            coherent(bus,line,meth,ns,tol)
% Syntax: [area,nmach_a,areal,nmach_al] = ...
%                           coherent(bus,line,meth,ns,tol)
%
% Purpose: Find the coherent machine groups using one of the
%          4 methods in PST.
%
% Input: bus - solved loadflow bus data
%        line - line data
%        meth - 1 for slow coherency, 2 for tight coherency
%               3 for Zaborszky, 4 for weak link, 
%               5 for Eliasson
%
% Output: area - matrix of coherent groups of machines
%         nmach_a - vector of number of machines in each 
%                   coherent area
%         areal - matrix of loose coherent groups
%         nmach_al - vector of number of machines in each
%                    loose coherent area
%
% See Also: V_slow, coh_map, L_group, ex_group, zabors, 
%           weaklink

% Algorithm: 
%
% Calls:
%
% Call By:

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1992

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

disp('executing function COHERENT')

if meth == 1         % slow coherency grouping
    disp('generate -inv(M)*K matrix')
    MK = svm_em(bus,line,1);
    disp('extract first ns eigenvalues')
    [eig_s,V_s] = v_slow(MK,ns);
    % use eigensubspace V_s to find coherent groups
    [area,nmach_a] = l_group(V_s);
  elseif meth == 2   % tight coherency grouping
    disp('generate -inv(M)*K matrix')
    MK = svm_em(bus,line,1);
    disp('extract first ns eigenvalues')
    [eig_s,V_s] = v_slow(MK,ns);
    % generate coherency map
    [grouping,c_map] = coh_map(V_s,tol);
    mesh(c_map)
    disp('press any key to continue for a 2-D plot'),pause
    view(2)
    title('Coherency Map')
    % use coherency map to find tight areas
    [area,nmach_a,areal,nmach_al] = ex_group(grouping);
  elseif meth == 3   % Zaborszky weak coupling method
    disp('generate -inv(M)*K matrix')
    MK = svm_em(bus,line,1);
    [area,nmach_a] = zabors(MK,tol);
  elseif meth == 4   % weak line method
    disp('generate -inv(M)*K matrix')
    MK = svm_em(bus,line,1);
    [area,nmach_a,areal,nmach_al] = weaklink(MK,tol);
  elseif meth == 5   % Eliasson grouping
    disp('generate -inv(M)*K matrix')
    MK = svm_em(bus,line,1);
    disp('extract first ns eigenvalues')
    [eig_s,V_s] = v_slow(MK,ns);
    % use sign to find tight areas
    [area,nmach_a] = sign_coh(V_s);
  else
    disp('COHERENT: invalid choice of method')
end
