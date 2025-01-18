clc; close all; clear all;

A = [-0.313, 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
B = [0.232; 0.0203; 0];
Q = diag([1 1 30]);
R = 1;

K = lqr(A,B,Q,R)

% S aleg eu cum vreau diferita de 0 daca se poate, nu prea influenteaza de
% exemplu eye(3)

S = eye(3);


P = S;
dt = 0.0001;
t = 0:dt:30;
N = length(t);

for k=N:-1:2
    P = P + dt*(A'*P+P*A+Q-P*B*inv(R)*B'*P);
    K(k,:) = inv(R)*B'*P;
end

plot(K); grid;

K1 = lqr(A,B,Q,R)

% valorile de la K1 sunt exact cele ce le vedem in plotul de mai sus









