clc; close all; clear all;

A = [-0.313, 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
B = [0.232; 0.0203; 0];
Q = diag([1 1 30]);
R = 1;

K = lqr(A,B,Q,R)

C = [0 0 1];

% asta e fara efect integrator
N = -inv(C*inv(A-B*K)*B)

%% cand introduc efectul integrator calculez A extins si B extins (vezi lab la pagina 4)
clc;
C = [0 0 1];
Ae = [A zeros(3,1); -C 0];
Be = [B;0];
Q = diag([1 1 1 30]);
R = 1;

K = lqr(Ae,Be,Q,R)
Kz = K(4)
Kx = K(1:3)
% K o sa aiba 4 valori: 3 se vor folosi pentru Kx din simulink si a 4-a
% valoare va fi Kz din simulink

%% cerinta 3 cu impunere de poli
clc;
A = [-0.313, 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
B = [0.232; 0.0203; 0];
Q = diag([1 1 30]);
R = 1;

K = lqr(A,B,Q,R)

C = [0 0 1];

N = -inv(C*inv(A-B*K)*B);

poli = eig(A-B*K)
poli = real(poli)*5 + imag(poli)*i;
L = place(A',C',poli)
%%
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









