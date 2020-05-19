function Cg = makeCg(bus,gen)
% MAKECG   Builds the sparse generator connection matrix for other
% computations.
%
%     Cg = MAKEBDC(MPC)
%     Cg = MAKEBDC(BUS, GEN)
%
%  Returns the sparse bus to generator connection matrix. The entry (i,j)
%  is 1 if the ith bus houses the jth generator.
%

%% extract from MPC if necessary
if nargin < 2
    mpc  = bus;
    bus  = mpc.bus;
    gen  = mpc.gen;
end
%% constants
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of generators

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeBdc: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%% build the sparse matrix
Cg = sparse(gen(:,GEN_BUS), (1:ng)', ones(ng, 1), nb, ng); 

end

