function [order] = getorder(mpc)
%getorder function to obtain the order of buses in a matrix using the
%standard ordering.
%   Standard ordering: slack bus(es) - generator buses - load buses
%   -----------------------------------------------------------------------
%   INPUT
%   mpc: power system case struct
%   -----------------------------------------------------------------------
%   OUTPUT
%   order: order of the buses to reorder computed matrices

%% Function definition
gen = find(mpc.bus(:,2)==2);
load = find(mpc.bus(:,2)==1);
slack = find(mpc.bus(:,2)==3);
order = [slack;gen;load];
end

