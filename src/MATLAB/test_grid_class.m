% test Grid_class.m work 
% you can load any MATPOWER case.
%     for example  loadcase('case2737sop');
clc
res = [];
clear;

% res :  fast_full fast_loop brute_full cont C1_isl C2_isl final 


runcase = loadcase('case39');
grid = Grid_class(runcase,'case39');


grid = grid.N_1_analysis(); % runs N-1 analysis
grid = grid.N_2_analysis('fast'); % runs fast N-2 contingency analysis
return
res(3,1) = grid.t_fast;
res(3,2) = grid.t_fast - grid.t_brute;
res(3,5) = Sz.r(grid.C1_isl);
res(3,6) = Sz.r(grid.C2_isl);
%grid = grid.N_2_analysis('bruteforce'); % runs bruteforce N-2 contingency analysis
%res(3,3) = grid.t_brute;
%res(3,4) = Sz.r(grid.brute_cont);


runcase = loadcase('case2746wop');
grid = Grid_class(runcase,'poland_winter_peak');


grid.N_1_analysis(); % runs N-1 analysis
grid = grid.N_2_analysis('fast'); % runs fast N-2 contingency analysis
res(1,1) = grid.t_fast;
res(1,2) = grid.t_fast - grid.t_brute;
res(1,5) = Sz.r(grid.C1_isl);
res(1,6) = Sz.r(grid.C2_isl);
%grid = grid.N_2_analysis('bruteforce'); % runs bruteforce N-2 contingency analysis
%res(1,3) = grid.t_brute;
%res(1,4) = Sz.r(grid.brute_cont);

runcase = loadcase('case2383wp');
grid = Grid_class(runcase,'poland_winter_peak');


grid = grid.N_1_analysis(); % runs N-1 analysis
grid = grid.N_2_analysis('fast'); % runs fast N-2 contingency analysis
res(2,1) = grid.t_fast;
res(2,2) = grid.t_fast - grid.t_brute;
res(2,5) = Sz.r(grid.C1_isl);
res(2,6) = Sz.r(grid.C2_isl);
%grid = grid.N_2_analysis('bruteforce'); % runs bruteforce N-2 contingency analysis
%res(2,3) = grid.t_brute;
%res(2,4) = Sz.r(grid.brute_cont);



