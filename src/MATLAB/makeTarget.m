function mpc = makeTarget(mpc,altmpc)
%MAKETARGET Creates the target matpower case for continuation power flow.
%
%   [TARGETMPC] = MAKETARGET(MPC,ALTMPC)
%   MPC: matpower case representing the base case for cpf
%   ALTMPC: matpower case with the contingency
%   TARGETMPC: output matpower case representing the target case for cpf
%
%   Evaluates and alters the power injections (generation and load) at the
%   buses in order to carry out the continuation power flow.
%   
%   A line trip is not represented as BR_STATUS=0 in this case. Instead it
%   is represented as some altered injections at other buses. In order to
%   evaluate the altered injections, the contingency is applied and islands
%   and isolated nodes are identified. The loads and generations in the
%   isolated nodes are made zero and for islands, they are altered through
%   the redispatch function.
%   


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Make generation/load zero in isolated buses
[~,isolated] = find_islands(altmpc);

% Handle isolated buses
[~,busind] = ismember(isolated,altmpc.bus(:,BUS_I));
[~,genind] = ismember(isolated,altmpc.gen(:,GEN_BUS));
mpc.bus(busind,PD)=0;
mpc.bus(busind,QD)=0;
mpc.gen(genind,PG)=0;
mpc.gen(genind,VG)=0;
mpc.bus(busind,VM)=0;

%% Redispatch generation and demands in each island
mpc_array = extract_islands(altmpc);
for i=1:length(mpc_array)
    submpc = redispatch(mpc_array{1,i});
    
    % Replace Pd values
    [~,busind] = ismember(submpc.bus(:,BUS_I),altmpc.bus(:,BUS_I));
    mpc.bus(busind,PD) = submpc.bus(:,PD);
    mpc.bus(busind,QD) = submpc.bus(:,QD);
    
    % Replace Pg values
    subgen = submpc.gen(:,GEN_BUS);
    simgen = altmpc.gen(:,GEN_BUS);
    genind = ismember(simgen,intersect(subgen,simgen,'stable'));
    mpc.gen(genind,PG) = submpc.gen(:,PG);
end
end


