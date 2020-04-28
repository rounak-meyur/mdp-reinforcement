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
% This evaluates the cascade of trips for an outage in each branch in the
% network. The cascade of trips is stored as a matrix of cells with each
% row listing branchies which are tripped at successive iterations. The
% branches trips at each iteration are listed in each column. DC Power flow
% is used in this study.
trip = cell(nl,1);
for k=1:nl
    iter = 1;
    trip{k,iter}=k;

    % Run DC power flow in steady state
    flogic = zeros(nl,1);
    simmpc = rundcpf(mpc,mpopt,fname);

    while true
        % Simulate the trip event
        simmpc.branch(trip{k,iter},BR_STATUS) = 0;

        % Insert code to handle islands without slack bus
        simmpc = handleIslands(simmpc);

        % Run dc power flow
        simmpc = rundcpf(simmpc,mpopt,fname);
        pf = simmpc.branch(:,PF);
        flogic(:,iter) = ~(abs(pf)<=flim);
        
        % Check for out of limit lines
        if ismember(1,flogic(:,iter))
            trip{k,iter+1} = find(flogic(:,iter)==1);
            iter = iter+1;
        else
            break
        end

    end
end