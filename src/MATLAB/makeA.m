function [A] = makeA(mpc)
%makeA Function to compute the branch bus incidence matrix
%
%   INPUT
%   mpc: power system case struct
%   -----------------------------------------------------------------------
%   OUTPUT
%   A: The branch-bus incidence matrix of the power system test case
%% Get bus incidence matrix
fbus = mpc.branch(:,1);
tbus = mpc.branch(:,2);
[~,fbusind] = ismember(fbus,mpc.bus(:,1));
[~,tbusind] = ismember(tbus,mpc.bus(:,1));
rows = [(1:nline)';(1:nline)'];
cols = [fbusind;tbusind];
entries = [ones(nline,1);-ones(nline,1)];
A = sparse(rows,cols,entries,nline,nbus);
end

