clc
clear

V_s = rand(16,5);

[nmach,narea] = size(V_s);
V = [zeros(1,narea+1); zeros(nmach,1) V_s];
mach_ord = [1:1:nmach]';
for i = 1:narea   % find reference machine for each area
  V = V(2:nmach-i+2,2:narea-i+2);
  [V_max,row_ind] = max(abs(V)); % find max in each column
                                 % for reference machine
  [V_max,col_ind] = max(V_max);  % find overall max
  % rearrange the machine ordering
  new_ref = row_ind(col_ind);
  temp = mach_ord(i);
  mach_ord(i) = mach_ord(new_ref+i-1);
  mach_ord(new_ref+i-1) = temp;
  % permute the V matrix
  temp = V(1,:);          % 1st do column permutation
  V(1,:) = V(new_ref,:); 
  V(new_ref,:) = temp;
  temp = V(:,1);          % next do row permutation
  V(:,1) = V(:,col_ind);
  V(:,col_ind) = temp;
  % perform Gaussian elimination
  v = V(:,1)/V(1,1);
  X = v*V(1,:);
  V = V - X;
end