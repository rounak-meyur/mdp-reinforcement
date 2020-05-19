function mpc = redispatch(mpc)
%redispatch Redispatches generation based on governor action and sheds load
%based on simple load shedding module.
%
%
%   -----------------------------------------------------------------------
%   Inputs:
%   mpc: input power system case
%   -----------------------------------------------------------------------
%   Outputs:
%   mpc: output power system case
%   -----------------------------------------------------------------------

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[XL, RA, XD, XDP, XDPP, TDP, TDPP, XQ, XQP, XQPP, TQP, TQPP,...
    H, D0, D1, S1, S2, RGOV, T1, T2, T3, T4, T5, F] = idx_gendyn;

%% Convert to internal ordering
intmpc = ext2int(mpc);
nb = size(intmpc.bus,1);
ng = size(intmpc.gen,1);

%% Compute the generation-load mismatch in the island
Cg = makeCg(intmpc);                    % Bus-generator connection matrix
LOSS = real(sum(get_losses(intmpc)));   % Total losses
LOAD = sum(intmpc.bus(:,PD));           % Total load in island
GEN = sum(Cg*intmpc.gen(:,PG));         % Total generation in island
MISMATCH = LOAD-GEN+LOSS;               % Total load-generation mismatch

%% Compute frequency deviation: positive implies drop, negative means increase
Rgov = intmpc.gen(:,RGOV);
Rinv = 1./Rgov;
Rinv(isinf(Rinv)) = 0;          % handle exception for uncontrolled generation
beta = sum(Rinv);               % stiffness
freqdev = MISMATCH/beta;

%% Redispatch generation or shed loads
if (-5.0<=freqdev)&&(freqdev<=0.8)
    % If frequency is in between 59.2Hz and 65.0 Hz, redispatch through
    % generator droop control
    delPg = freqdev*Rinv;
    delPd = zeros(nb,1);
    delQd = zeros(nb,1);
    disp("Generators redispatched...\n")
elseif freqdev<-5.0
    % If frequency is greater than 65.0 Hz, trip the generator over 
    % frequency relays
    delPg = intmpc.gen(:,PG);
    delPd = zeros(nb,1);
    delQd = zeros(nb,1);
    disp("Generators tripped by overfrequency relays.\n")
elseif (0.8<freqdev)&&(freqdev<=1.2)
    % If frequency is less than 59.2 Hz and more than 58.8 Hz, 10% of the 
    % load is removed
    delPg = zeros(ng,1);
    delPd = 0.1*intmpc.bus(:,PD);
    delQd = 0.1*intmpc.bus(:,QD);
    disp("10% of total load shed.")
elseif (1.2<freqdev)&&(freqdev<=2.0)
    % If frequency is less than 58.8 Hz and more than 58.0 Hz, 25% of the 
    % load is removed
    delPg = zeros(ng,1);
    delPd = 0.25*intmpc.bus(:,PD);
    delQd = 0.25*intmpc.bus(:,QD);
    disp("25% of total load shed.")
elseif (2.0<freqdev)&&(freqdev<=5.0)
    % If frequency is less than 58.0 Hz and more than 55.0 Hz, 45% of the 
    % load is removed
    delPg = zeros(ng,1);
    delPd = 0.45*intmpc.bus(:,PD);
    delQd = 0.45*intmpc.bus(:,QD);
    disp("45% of total load shed.")
else
    delPg = zeros(ng,1);
    delPd = intmpc.bus(:,PD);
    delQd = intmpc.bus(:,QD);
    disp("Frequency too low. All load shed.")
end

intmpc.bus(:,PD) = intmpc.bus(:,PD) - delPd;
intmpc.bus(:,QD) = intmpc.bus(:,QD) - delQd;
intmpc.gen(:,PG) = intmpc.gen(:,PG) - delPg;

%% Convert to original ordering
mpc = int2ext(intmpc);
end

