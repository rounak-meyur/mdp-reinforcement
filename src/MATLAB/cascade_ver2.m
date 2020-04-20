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
trip = cell(nl,1);
for k=1:nl
    iter = 1;
    trip{k,iter}=k;

    % Run AC power flow in steady state
    flogic = zeros(nl,1);
    simmpc = runpf(mpc,mpopt,fname);

    while true
        % Simulate the trip event
        simmpc.branch(trip{k,iter},BR_STATUS) = 0;

        % Insert code to handle islands without slack bus
        simmpc = handleIslands(simmpc);

        % Run ac power flow
        simmpc = runpf(simmpc,mpopt,fname);
        pf = simmpc.branch(:,PF);
        qf = simmpc.branch(:,QF);
        sf = abs(pf + 1j*qf);
        flogic(:,iter) = ~(abs(sf)<=flim);
        
        % Check for out of limit lines
        if ismember(1,flogic(:,iter))
            trip{k,iter+1} = find(flogic(:,iter)==1);
            iter = iter+1;
        else
            break
        end

    end
end