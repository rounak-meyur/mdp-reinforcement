clc
clear

BR_STATUS = 11;
RATE_A = 6;
PF = 14;
mpc = loadcase('case39');

%% Option override
mpopt = mpoption;
mpopt = mpoption(mpopt, 'verbose', 0);
fname = 'power-flow.txt';

%% Initialize simulation
% Run DC power flow in steady state
simmpc = rundcpf(mpc,mpopt,fname);
return
flim = simmpc.branch(:,RATE_A);
f = simmpc.branch(:,PF);
flogic = abs(f)<=flim;

%% Event simulation (line trip with no transient)
% Simulate the trip event
nl = size(mpc.branch,1);
for k = 1:14
    simmpc = rundcpf(mpc,mpopt,fname);
    simmpc.branch(k,BR_STATUS) = 0;

    % Insert code to handle islands without slack bus
    [groups,isolated] = find_islands(simmpc);
    mpc_array = extract_islands(simmpc,groups);
    size(mpc_array)
end
return
% Run dc power flow
simmpc = rundcpf(simmpc);

% Check power flow
flim = simmpc.branch(:,RATE_A);
f = simmpc.branch(:,PF);
flogic = abs(f)<=flim;

%% PTDF based computations
% Check flow from PTDF
S = makePTDF(simmpc);
s = makeSbus(simmpc.baseMVA,simmpc.bus,simmpc.gen);
p = real(s);
fc = simmpc.baseMVA*S*p;