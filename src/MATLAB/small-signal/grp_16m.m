% grp_16m.m
% m.file to find coherent machine groups in the 16
%     machine, 68 bus system in LOTDYS report 
%     using the Matlab Power System Toolbox
%
% generate the inv(M)*K matrix using linearization
%svm_16mk
%save svm_16mk MK

load svm_16mk

% extract first five eigenvalues
[eig_s,V_s] = V_slow(MK,5);

[area,nmach_a,L] = L_group(V_s);

%tol = 0.9
[grouping,c_map] = coh_map(V_s,tol);
[nx,areax,nmach_x,areax_l] = ex_group(grouping);
m = [30 60];  % orientation for 3D plot
