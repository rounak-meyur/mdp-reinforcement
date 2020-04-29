function [grouping,c_map] = coh_map(V_s,tol)
% Syntax: [grouping,c_map] = coh_map(V_s,tol)
%
% Purpose: Set up the coherency map 
%
% Input: V_s - eigenvector matrix of slow eigenvalue
%        tol - tolerance for coherency
% Output: grouping - a matrix of the cosine of the angles  
%                    between the machines - the tolerance
%         c_map - the grouping matrix with negative elements
%                 set to zero
%
% See Also: V_slow, L_group

% Algorithm: Compute the cosine of the angle between the 
%            machine mode shapes for the slow eigenvalues.
%            A number close to 1 indicates strong coherency.
%            The matrix grouping is best display with the 
%            format +.
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

[n narea] = size(V_s);
% normalize the V_s columns
for i = 1:narea
  V_s(:,i) = V_s(:,i)/norm(V_s(:,i));
end

% normalize row of V_s
for i = 1:n
 V_s(i,:) = V_s(i,:)/norm(V_s(i,:));
end

% compute coherency map
c_map = V_s*V_s'; 
grouping = c_map - tol*ones(n,n);
c_map = max(grouping,zeros(n,n));

