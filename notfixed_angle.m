function f = notfixed_angle(x)
load('G.mat');
load('mean_deep.mat');
n = x(1);
q_eq = x(2:4);
q_eq = reshape(q_eq,3,1);
GQ = G(:,3*n-2:3*n) * q_eq;

f = norm(M2 - GQ,2);