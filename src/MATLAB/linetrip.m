function [base,target] = linetrip(mpc,brnchind)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% information of tripped lines
mpc.branch(brnchind,BR_STATUS)=0;
fbus = mpc.branch(brnchind,F_BUS);
tbus = mpc.branch(brnchind,T_BUS);
pf = mpc.branch(brnchind,PF);
qf = mpc.branch(brnchind,QF);
pt = mpc.branch(brnchind,PT);
qt = mpc.branch(brnchind,QT);

%% Handle islands to reassign bus types
mpc = handleIslands(mpc);

%% base matpower case: replace the tripped line(s) with equivalent line injections
base = mpc;
[~,fbusind] = ismember(fbus,base.bus(:,BUS_I));
[~,tbusind] = ismember(tbus,base.bus(:,BUS_I));

for i=1:length(brnchind)
    base.bus(fbusind(i),PD) = base.bus(fbusind(i),PD)+pf(i);
    base.bus(fbusind(i),QD) = base.bus(fbusind(i),QD)+qf(i);
    base.bus(tbusind(i),PD) = base.bus(tbusind(i),PD)+pt(i);
    base.bus(tbusind(i),QD) = base.bus(tbusind(i),PD)+qt(i);
end
%% target matpower case: redispatch generation and shed loads
target = mpc;
mpc_array = extract_islands(mpc);
for i=1:length(mpc_array)
    submpc = redispatch(mpc_array{1,i});
    
    % Replace Pd values
    [~,busind] = ismember(submpc.bus(:,BUS_I),mpc.bus(:,BUS_I));
    target.bus(busind,PD) = submpc.bus(:,PD);
    target.bus(busind,QD) = submpc.bus(:,QD);
    
    % Replace Pg values
    subgen = submpc.gen(:,GEN_BUS);
    simgen = mpc.gen(:,GEN_BUS);
    genind = ismember(simgen,subgen);
    target.gen(genind,PG) = submpc.gen(:,PG);
end


end

