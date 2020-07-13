clc
clear

c = ones(5,1);
C = ones(5);
z = zeros(5,1);
W = eye(5)-0.2*C;
Q = sqrtm(W);
X = sdpvar(5,5,'full');
F = sdpvar(5);

M = [X Q;Q F+0.2*C];

objective = trace(X);
constraints = [M>=0, F*c==z,X>0,F>0];
opt = sdpsettings('debug','1','solver','sedumi','dualize','1');

sol = optimize(constraints,objective,opt);

if sol.problem==0
    Fval = value(F);
    Xval = value(X);
else
    fprintf("No solution");
end

return