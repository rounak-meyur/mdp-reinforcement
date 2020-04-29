function [eig_s,V_s] = V_slow(MK,n_s)
% Syntax: [eig_s,V_s] = V_slow(MK,n_s)
%
% Purpose: Set up the slow eigensubspace for coherent
%          machine grouping
%
% Input: MK - matrix inv(M)*K
%        n_s - number of slow modes
%
% Output: eig_s - vector of slow eigenvalues
%         V_s - eigenvectors corresponding to the slow
%               eigenvalues
%
% See Also: L_group

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     May 1991

[eig_vec,lambda] =eig(MK);
[lambda,k] = sort(abs(diag(lambda)));
eig_s = lambda(1:n_s);
V_s = eig_vec(:,k(1:n_s));
