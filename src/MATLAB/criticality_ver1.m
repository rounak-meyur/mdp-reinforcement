clc
clear

%% Load the case
%workdir = pwd;
%mpc = ACTIVSg2000;
mpc = loadcase('case39');

%% Create a test scenario
% Sample power injections
p = sampleinj(mpc);

% Rating of lines
flim = 1.5*mpc.branch(:,6)/mpc.baseMVA;


% Check random draws
%check_invalid = (pgmin>pg) | (pg>pgmax);

% PTDF matrix
k = 27;
mpc.branch(k,11) = 0;
% mpc.bus(33,2)=3;
mpc.branch(:,9) = ones(46,1);
mpc = ext2int(mpc);
S = makePTDF(mpc);
% mpc = int2ext(mpc);
% flim(k) = [];
% f = S*p;
% checkflow = (abs(f) <= flim);
% find(checkflow==0)