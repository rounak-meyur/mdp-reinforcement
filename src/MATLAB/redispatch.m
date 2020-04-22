function [mpc] = redispatch(mpc)
%redispatch Optimization framework to shed load based on the DC OPF
%formulation.
%
%   The function takes as input the current generations and loads at all
%   buses and solves the optimization problem. The output is the new set of
%   generations and loads at all buses.
%
%   -----------------------------------------------------------------------
%   Inputs:
%   mpc: input power system case
%   -----------------------------------------------------------------------
%   Outputs:
%   mpc: output power system case
%   -----------------------------------------------------------------------

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Check if all buses are disconnected
btype = mpc.bus(:,BUS_TYPE);
if min(btype)==4 && max(btype==4)
    mpc.bus(:,PD)=zeros(length(btype),1);
    mpc.bus(:,QD)=zeros(length(btype),1);
else
    %% Convert to internal ordering
    intmpc = ext2int(mpc);

    %% System constants
    nb = size(intmpc.bus,1);

    % Branch data
    fmax = intmpc.branch(:,RATE_A)/intmpc.baseMVA;

    % Generator data
    ng = size(intmpc.gen,1);
    gmax = intmpc.gen(:,PMAX)/intmpc.baseMVA;
    gmin = intmpc.gen(:,PMIN)/intmpc.baseMVA;
    pg = intmpc.gen(:,PG)/intmpc.baseMVA;

    % Generator to Bus Connection Matrix
    BG = sparse(intmpc.gen(:,GEN_BUS),(1:ng)',ones(ng,1),nb,ng);

    % Load Data
    dmax = intmpc.bus(:,PD)/intmpc.baseMVA;
    dmin = zeros(nb,1);

    % PTDF matrix
    S = makePTDF(intmpc);

    %% Optimization for load shedding
    g = sdpvar(ng,1);
    d = sdpvar(nb,1);

    c1 = ones(nb,1);
    c2 = ones(ng,1);

    constraints = [gmin<=g<=gmax, dmin<=d<=dmax,...
        c1'*((BG*g)-d)==0, -fmax<=S*((BG*g)-d)<=fmax];
    objective = c1'*(dmax-d)+c2'*abs(pg-g);

    optimize(constraints,objective);
    g_solved = value(g);
    d_solved = value(d);

    %% Update the power system case
    ratio = intmpc.bus(:,QD)./intmpc.bus(:,PD);
    ratio(isnan(ratio))=0;
    intmpc.bus(:,PD) = d_solved * intmpc.baseMVA;
    intmpc.bus(:,QD) = intmpc.bus(:,PD) .* ratio;
    intmpc.gen(:,PG) = g_solved * intmpc.baseMVA;

    mpc = int2ext(intmpc);
end
end

