% grp_ipp.m
% script file to find coherent machine groups in the 329
%     machine IPP system 
%     using the Matlab Power System Toolbox
%
% load eigenvector matrix of first 15 modes
ipp15

%[area,nmach_a,L] = L_group(V_s);

while(0)
clock1 = clock;
tol = 0.95 % cos(18.2deg) = 0.95
[grouping] = coh_map(L,tol);
[nx_9,areax_9,nmach_x_9] = ex_group(grouping);
%m = [30 60];
t9 = etime(clock,clock1);
save ipp_9
end %while(0)

while(0)
clock1 = clock;
tol = 0.925 % cos(22.3deg) = 0.925
[grouping] = coh_map(L,tol);
[nx_925,areax_925,nmach_x_925] = ex_group(grouping);
%m = [30 60];
t925 = etime(clock,clock1);
end

clock1 = clock;
tol = 0.985 % cos(10deg) = 0.9848
% eliminate the first column of L as an experiment
[nrow ncol] = size(L);
L = L(:,2:ncol);
[grouping] = coh_map(L,tol);
[nx_985,areax_985,nmach_x_985,loose_985] = ex_group(grouping);
%m = [30 60];
t985 = etime(clock,clock1);
save ipp_985

