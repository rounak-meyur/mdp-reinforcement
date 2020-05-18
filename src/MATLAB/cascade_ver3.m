clc
clear

%% Setup MATLAB directories
cd ..
pardir = pwd;
prgdir = strcat(pardir,'\progress\');
figdir = strcat(pardir,'\figures\');
casdir = strcat(pardir,'\case-data\');
cd .\MATLAB;

%% Constants/Input information
BR_STATUS = 11;
RATE_A = 6;
PF = 14;
QF = 15;

%% Option override
mpopt = mpoption('verbose',0,'out.all',0);
suffix = datestr(now,'mm-dd-yy-HH-MM-SS-FFF');
fname = strcat(prgdir,'power-flow',suffix,'.txt');

%% Initialize simulation
% Load case and 
mpc = loadcase('case39');
mpc = case39mac(mpc);
nl = size(mpc.branch,1);
flim = mpc.branch(:,RATE_A);


%% Event simulation (line trip with no transient)
%   Steady state power flow: Newton Raphson load flow model
%   Initiating event as the first set of lines
%   List of steps in the cascading model
%   1. Trip out of limit transmission lines/transformers
%   2. Redispatch generators
%   3. Loadshedding
%   4. Continuation power flow to identify voltage instability
%   


trip_k = 27;

% Step 0: Run AC power flow in steady state
flogic = zeros(nl,1);
simmpc = runpf(mpc,mpopt,fname);

% Start of the cascading scenario
% Step 1: Simulate the trip event
simmpc.branch(trip_k,BR_STATUS) = 0;

% Step 2: Find islands and compute load-generation imbalance
[groups,isolated] = find_islands(simmpc);
% Isolated nodes are to be disconnected
% For each island compute the load-generation imbalance
% Compute the frequency deviation based on the droop characteristics
% 
Cg = makeCg(simmpc);

return
for j=1:length(groups)
    [~,busind] = ismember(groups{1,j},simmpc.bus(:,1)); % index of buses
    bustype = simmpc.bus(busind,2);                     % type of bus
    
    % Compute total load, generation and mismatch
    LOAD = sum(simmpc.bus(busind,3));
end


return
%% Test setup
% Run ac power flow after redispatch
simmpc_target = runpf(simmpc_mod,mpopt,fname);
% pf = simmpc.branch(:,PF);
% qf = simmpc.branch(:,QF);
% sf = abs(pf + 1j*qf);
% flogic(:,iter) = ~(abs(sf)<=flim);

%% Continuation power flow
simmpc_init.bus(:,2) = simmpc_mod.bus(:,2);
results = runcpf(simmpc_init,simmpc_target,mpopt,fname);

vmag = abs(results.cpf.V);
vang = angle(results.cpf.V);

subplot(1,2,1)
plot(vmag(30:39,:)')
subplot(1,2,2)
plot(vang(30:39,:)')

