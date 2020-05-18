function [mpc] = handleIslands(mpc)
%HANDLEISLANDS Handles islands in a network by setting slack buses.
%
%   [OUTMPC] = HANDLEISLANDS(MPC)
%   Scans through the network to identify isolated nodes and islands in the
%   network. Returns a Matpower struct with the following considerations.
%   
%   a.  isolated nodes are considered to be disconnected
%   b.  if an island has no PV/slack bus, all nodes are considered to be
%   disconnected.
%   c.  if an island has no slack bus, but has at least one pv bus,
%   identify a suitable slack bus among them.
%   


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Handle islands and isolated buses
[groups,isolated] = find_islands(mpc);

% Handle islands
for j=1:length(groups)
    [~,busind] = ismember(groups{1,j},mpc.bus(:,BUS_I));    % index of buses
    bustype = mpc.bus(busind,BUS_TYPE);                     % type of bus
    
    if ismember(REF,bustype)==0
        % If no slack bus in the group
        if ismember(PV,bustype)==0
            % If no PV bus, disconnect the buses
            mpc.bus(busind,BUS_TYPE)=NONE;
        else
            % If at least one PV bus
            pvind = busind(bustype==PV);
            pvbus = mpc.bus(pvind,BUS_I);
            % Compute total generation capacity on each PV bus
            pmax = zeros(length(pvbus),1);
            for i=1:length(pvbus)
                % index of generators on pv bus
                genind = mpc.gen(:,GEN_BUS)==pvbus(i);
                pmax(i) = sum(mpc.gen(genind,PMAX));
            end
            busind = pvind(pmax==max(pmax));
            mpc.bus(busind,BUS_TYPE)=REF;
        end
    end
end

% Handle isolated buses
[~,busind] = ismember(isolated,mpc.bus(:,BUS_I));
mpc.bus(busind,BUS_TYPE)=4;

%% Redispatch generation and demands in each island
mpc_array = extract_islands(mpc);
for i=1:length(mpc_array)
    submpc = optredispatch(mpc_array{1,i});
    
    % Replace Pd values
    [~,busind] = ismember(submpc.bus(:,BUS_I),mpc.bus(:,BUS_I));
    mpc.bus(busind,PD) = submpc.bus(:,PD);
    mpc.bus(busind,QD) = submpc.bus(:,QD);
    
    % Replace Pg values
    subgen = submpc.gen(:,GEN_BUS);
    simgen = mpc.gen(:,GEN_BUS);
    genind = ismember(simgen,intersect(subgen,simgen,'stable'));
    mpc.gen(genind,PG) = submpc.gen(:,PG);
end
end

