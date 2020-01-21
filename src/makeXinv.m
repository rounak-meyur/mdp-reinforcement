function [Xinv, b] = makeXinv(baseMVA, bus, branch)
%MAKEBDC   Builds the inverse X matrix and vector of susceptances.
%   [Xinv, b] = MAKEX(MPC)
%   [Xinv, b] = MAKEX(BASEMVA, BUS, BRANCH)
%
%   Returns the inverse X matrix and susceptance vectors. 

%   Xinv matrix is a digonal matrix with entries as the brance susceptances.
%   b is a vector of the branch susceptances.

%   Does appropriate conversions to p.u.
%   Bus numbers must be consecutive beginning at 1 (i.e. internal ordering).
%
%   Example:
%       [Xinv, b] = makeX(baseMVA, bus, branch);

%% extract from MPC if necessary
if nargin < 3
    mpc     = baseMVA;
    baseMVA = mpc.baseMVA;
    bus     = mpc.bus;
    branch  = mpc.branch;
end

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeBdc: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

%% compute the vector of branch susceptances.
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;

%% build X
Xinv = diag(b);
