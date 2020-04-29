% grp_48m.m
% script file to find coherent machine groups in the 48
%     machine NPCC system 
%     using the Matlab Power System Toolbox
%

l
V_s = L;

[area,nmach_a,L_x] = L_group(V_s);

%tol = 0.9
[nrow ncol] = size(V_s);
[grouping,c_map] = coh_map(V_s,tol);
t1 = clock;
[nx,areax,nmach_x,area_l] = ex_group(grouping);
etime(clock,t1)
m = [30 60];  % orientation for 3D plot