clc
clear

%% Setup MATLAB directories
cd ..
pardir = pwd;
prgdir = strcat(pardir,'\progress\');
figdir = strcat(pardir,'\figures\');
casdir = strcat(pardir,'\case-data\');
cd .\MATLAB;



%% Option override
mpopt = mpoption('verbose',0,'out.all',0);
suffix = datestr(now,'mm-dd-yy-HH-MM-SS-FFF');
fname = strcat(prgdir,'power-flow',suffix,'.txt');

%% Initialize simulation
% Load case and 
mpc = loadcase('case39');
mpc = case39mac(mpc);
nl = size(mpc.branch,1);
flim = mpc.branch(:,6);


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
altmpc = simmpc;
altmpc.branch(trip_k,11) = 0;
simmpc_mod = makeTarget(simmpc,altmpc);


%% Continuation power flow
results = runcpf(simmpc,simmpc_mod,mpopt,fname);
vmag = abs(results.cpf.V);
plot(vmag')
grid on
check_vmag = diff(vmag');

