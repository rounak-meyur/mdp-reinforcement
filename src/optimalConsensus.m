clc
clear

c = ones(5,1);
C = ones(5);
z = zeros(5,1);
I = eye(5);
W = I-0.2*C;
Q = sqrtm(W);
X = sdpvar(5,5,'full');
F = sdpvar(5);

M = [X I;I F];

objective = trace(W*X);
constraints = [M>=0, X>=I, F>=I];

opt = sdpsettings('debug','1','solver','sedumi','dualize','1');

sol = optimize(constraints,objective,opt);

if sol.problem==0
    Fval = value(F);
    Xval = value(X);
else
    fprintf("No solution");
end

return