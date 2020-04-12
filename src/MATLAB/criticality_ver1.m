clc
clear

% Load the case
workdir = pwd;
mpc = case_ACTIVSg2000;
%mpc = loadcase('case39');

% Sample power injections
p = sampleinj(mpc);

% Rating of lines
flim = 1.5*mpc.branch(:,6)/mpc.baseMVA;


% Check random draws
check_invalid = (pgmin>pg) | (pg>pgmax);

% PTDF matrix
k = 6;
mpc.branch(k,11) = 0;
mpc = ext2int(mpc);
S = makePTDF(mpc);
mpc = int2ext(mpc);
flim(k) = [];
f = S*pinj;
checkflow = (abs(f) <= flim);
find(checkflow==0)