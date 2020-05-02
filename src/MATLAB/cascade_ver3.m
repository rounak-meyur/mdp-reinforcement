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
nl = size(mpc.branch,1);
flim = mpc.branch(:,RATE_A);


%% Event simulation (line trip with no transient)
trip_k = 20;

% Run AC power flow in steady state
flogic = zeros(nl,1);
simmpc_init = runpf(mpc,mpopt,fname);


% Simulate the trip event
simmpc_init.branch(trip_k,BR_STATUS) = 0;

% Handle islands without slack bus and
% redispatch load to maintain a healthy power system
simmpc_mod = handleIslands(simmpc_init);

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

