function [p] = sampleinj(mpc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Number of entries
nb = size(mpc.bus,1);

% Range of power demands
pdmin = zeros(nb,1);
pdmax = mpc.bus(:,3)/mpc.baseMVA;

% Get the generators
busnum = mpc.bus(:,1);
genbus = mpc.gen(:,1);
ng = size(mpc.gen,1);

% Range of generations
pgmin = zeros(nb,1);
pgmax = zeros(nb,1);
for i=1:ng
    busind = find(busnum==genbus(i));
    pgmin(busind) = pgmin(busind) + mpc.gen(i,10);
    pgmax(busind) = pgmax(busind) + mpc.gen(i,9);
end
pgmin = pgmin/mpc.baseMVA;
pgmax = pgmax/mpc.baseMVA;

% Draw random demands
pd = (pdmax-pdmin) .* rand(nb,1) + pdmin;
K = sum(pd-pgmin)/sum(pgmax-pgmin);
pg = pgmin + K*(pgmax-pgmin);
p = pg - pd;
end