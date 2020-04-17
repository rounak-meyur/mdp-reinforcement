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
TYPE = 2;
BUS = 1;
PF = 14;

%% Option override
mpopt = mpoption('verbose',0,'out.all',0);
suffix = datestr(now,'mm-dd-yy-HH-MM-SS-FFF');
fname = strcat(prgdir,'power-flow',suffix,'.txt');

%% Initialize simulation
% Load case and run DC power flow in steady state
mpc = loadcase('case39');
simmpc = rundcpf(mpc,mpopt,fname);

% Check power flow
flim = simmpc.branch(:,RATE_A);
f = simmpc.branch(:,PF);
flogic = abs(f)<=flim;

%% Event simulation (line trip with no transient)
% Simulate the trip event
nl = size(mpc.branch,1);
for k = 27:27
    simmpc = rundcpf(mpc,mpopt,fname);
    simmpc.branch(k,BR_STATUS) = 0;

    % Insert code to handle islands without slack bus
    [groups,isolated] = find_islands(simmpc);
    mpc_array = extract_islands(simmpc,groups);
    
    for j=1:size(groups,2)
        submpc = mpc_array{1,j};
        busnum = groups{1,j};
        btype = submpc.bus(:,TYPE);
        if ismember(3,btype)==0
            if ismember(2,btype)==0
                % If no PV bus, disconnect the buses
                [~,busid] = ismember(busnum,simmpc.bus(:,BUS));
                simmpc.bus(busid,TYPE)=4;
            else
                [~,subbusid]=ismember(2,btype);
                [~,busid]=ismember(busnum(subbusid),simmpc.bus(:,BUS));
                simmpc.bus(busid,TYPE)=3;
            end
        end
    end
    
    % Handle isolated buses
    [~,busid] = ismember(isolated,simmpc.bus(:,BUS));
    simmpc.bus(busid,TYPE)=4;
    
    % Run dc power flow
    simmpc = rundcpf(simmpc,mpopt,fname);
end
return

% Check power flow
f = simmpc.branch(:,PF);
flogic = abs(f)<=flim;

