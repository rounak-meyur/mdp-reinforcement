clc
clear

% Load the MATPOWER case using designated casefilenames
casefilename = 'case39';
mpc = loadcase(casefilename);
mpopt = mpoption('model', 'DC');
fname = 'dcpf.txt';
[results, success] = rundcpf(mpc, mpopt, fname);
% define_constants;

% LODF calculation for multiple contingencies
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(mpc);
[Ainc] = makeIncidence(mpc);
[Xinv, b] = makeXinv(mpc);

% Reduced matrices
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
Bred = Bbus([pv; pq], [pv; pq]);
Ared = Ainc(:, [pv; pq]);

% Get power flow results
f = abs(results.branch(:,14)/results.baseMVA);
f_lim = results.branch(:,6)/results.baseMVA;

% Line outage result
alpha = 4;          % first outaged line
beta = 10;           % second outaged line

% Line outage distribution factors for outage of line alpha
a_alpha = Ared(alpha,:);
dist_alpha = (b .* (Ared / Bred) * a_alpha') / ...
    (1 - b(alpha) * (a_alpha / Bred) * a_alpha');

% Line outage distribution factors for outage of line beta
a_beta = Ared(beta,:);
dist_beta = (b .* (Ared / Bred) * a_beta') / ...
    (1 - b(beta) * (a_beta / Bred) * a_beta');

% Get the double line outage terms
tau_alphabeta = (1+(dist_alpha(beta)*f(beta)/f(alpha)))/...
    (1-(dist_alpha(beta)*dist_beta(alpha)));
tau_betaalpha = (1+(dist_beta(alpha)*f(alpha)/f(beta)))/...
    (1-(dist_beta(alpha)*dist_alpha(beta)));

epsilon_alpha = (dist_alpha ./ (f_lim - f)) * f(alpha);
epsilon_beta = (dist_beta ./ (f_lim - f)) * f(beta);

constraints = epsilon_alpha*tau_alphabeta + epsilon_beta*tau_betaalpha;

% Clear screen
clc