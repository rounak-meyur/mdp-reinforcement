function [mpc] = loadshed(mpc)
%Loadshed Optimization framework to shed load based on the DC OPF
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

%% System inputs
nbus = length(mpc.bus);

% Branch data
fmax = mpc.branch(:,6)/mpc.baseMVA;

% Generator data
[~,genind] = ismember(mpc.gen(:,1),mpc.bus(:,1));
pg_current = zeros(nbus,1);
pg_current(genind,1) = mpc.gen(:,2)/mpc.baseMVA;

% Load Data
pd_current = mpc.bus(:,3)/mpc.baseMVA;

% PTDF matrix
S = makePTDF(mpc);
pg_current(all(isnan(S),1)) = 0;
mpc.bus(all(isnan(S),1),2) = 3;

%% Optimization for load shedding
pg = sdpvar(nbus,1);
pd = sdpvar(nbus,1);
c = ones(nbus,1);
z = zeros(nbus,1);

constraints = [z<=pg<=pg_current, z<=pd<=pd_current,...
    c'*(pg-pd)==0, -fmax<=S*(pg-pd)<=fmax];
objective = c'*(pd_current-pd);

optimize(constraints,objective);
pg_solved = value(pg);
pd_solved = value(pd);

%% Update the power system case
mpc.bus(:,3) = pd_solved * mpc.baseMVA;
mpc.gen(:,2) = pg_solved(genind,1) * mpc.baseMVA;
end

