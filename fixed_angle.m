function f = fixed_angle(x)
load('G.mat');
load('LocMat.mat');
load('mean_superficial.mat');
n = x(1);
q = x(2);
e_q = LocMat(:,n) ./ norm(LocMat(:,n),2);
GQ = G(:,3*n-2:3*n) * (q*e_q);

f = norm(M - GQ,2);