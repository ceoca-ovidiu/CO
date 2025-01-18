%% Rotary Inverted Pendulum
clc; close all; clear all;

m = 0.027;       % Mass of pendulum [kg]
g = 9.81;        % Gravitation acceleration [m/s^2]

Kt = 0.0531;     % Torque constant DC motor
Rm = 11.7356;    % Armature resistance

d1 = Kt/Rm;
d2 = (Kt^2)/Rm;

L2 = 0.328;      % Length of pendulum [m]
J2 = 0.0046617;  % Inertial moment of pendulum [kgm^2]
C2 = 0.0017;     % Friction coefficient of pendulum [Nm-s]

L1 = 0.205;        % Length of arm [m]
J1 = 0.0019;     % Inertial moment of arm [kgm^2]
C1 = 0.025;      % Friction coefficient of arm [Nm-s]

lambda = (J1 + m*(L1^2))*(J2 + m*(L2^2)) - (m*L1*L2)^2;

A21 = m*g*L2*(J1+m*(L1^2));
A22 = -C2*(J1+m*(L1^2));
A23 = 0;
A24 = -(m*L1*L2)*(C1+d2);

A41 = (m^2)*g*L1*(L2^2);
A42 = -C2*(m*L1*L2);
A43 = 0;
A44 = (J2+m*(L2^2))*(C1+d2);

B21 = d1*(m*L1*L2);
B41 = d1*(J2+m*(L2^2));

A = [0 lambda 0 0;
     A21 A22 A23 A24;
     0 0 0 lambda;
     A41 A42 A43 A44]./lambda

B = [0; B21; 0; B41]./lambda

if rank(ctrb(A,B))~=4
    error('Sistemul nu e controlabil');
end

% C = [1 0 0 0;
%      0 0 1 0];

%% Analiza sistemului in bucla deschisa
clc; close all;

eigValues = eig(A)

stepFinalValue = 1;
initialConditions = [1;1;1;1];
simulinkOut = sim('simulare_bucla_deschisa.slx');

time = simulinkOut.openLoopSimulinkOut.Time;
figure('Position', [550, 550, 900, 600]);
subplot(4,1,1)
plot(time, simulinkOut.openLoopSimulinkOut.Data(:,1)); grid; xlabel('Timp'); title('X1');
subplot(4,1,2)
plot(time, simulinkOut.openLoopSimulinkOut.Data(:,2)); grid; xlabel('Timp'); title('X2');
subplot(4,1,3)
plot(time, simulinkOut.openLoopSimulinkOut.Data(:,3)); grid; xlabel('Timp'); title('X3');
subplot(4,1,4)
plot(time, simulinkOut.openLoopSimulinkOut.Data(:,4)); grid; xlabel('Timp'); title('X4');

%% LQR pentru probleme de stabilizare (K FIX)
clc; close all;
Q = eye(4);
R = 1;
K_fix = lqr(A,B,Q,R)

simulinkOut = sim('LQR_probleme_de_stabilizare.slx')

time = simulinkOut.simoutResponse.Time;
figure('Position', [550, 550, 900, 600]);

subplot(4,1,1)
plot(time, simulinkOut.simoutResponse.Data(:,1)); grid; title('X1'); xlabel('Timp');
subplot(4,1,2)
plot(time, simulinkOut.simoutResponse.Data(:,2)); grid; title('X2'); xlabel('Timp');
subplot(4,1,3)
plot(time, simulinkOut.simoutResponse.Data(:,3)); grid; title('X3'); xlabel('Timp');
subplot(4,1,4)
plot(time, simulinkOut.simoutResponse.Data(:,4)); grid; title('X4'); xlabel('Timp');

figure('Position', [550, 550, 900, 600]);
plot(time, simulinkOut.simoutCommand.Data(:,1)); grid; title('Comandă'); xlabel('Timp');
%% LQR pentru probleme de stabilizare (K VARIABIL)
clc; close all;
Q = eye(4);
R = 1;
S = eye(4);
P = S;
dt = 0.0001;
t = 0:dt:30;
N = length(t);

for k=N:-1:2
    P = P + dt*(A'*P+P*A+Q-P*B*inv(R)*B'*P);
    K(k,:) = inv(R)*B'*P;
end

figure('Position', [550, 550, 900, 600]);
plot(K); grid;


x=[0 0 20 0]';
u=0;
for k=1:N-1
    u(k)=-K(k,:)*x(:,k);
    K1=sis_liniar(x(:,k),u(k));
    K2=sis_liniar(x(:,k)+0.5*dt*K1,u(k));
    K3=sis_liniar(x(:,k)+0.5*dt*K2,u(k));
    K4=sis_liniar(x(:,k)+dt*K3,u(k));
    
    x(:,k+1)=x(:,k)+1/6*dt*(K1+2*K2+2*K3+K4);   
end
figure('Position', [550, 550, 900, 600]);
subplot(4,1,1)
plot(x(1,:),'LineWidth',1)
title('x1'); xlabel('Timp')
grid
subplot(4,1,2)
plot(x(2,:),'LineWidth',1)
title('x2'); xlabel('Timp')
grid
subplot(4,1,3)
plot(x(3,:),'LineWidth',1)
title('x3'); xlabel('Timp')
grid
subplot(4,1,4)
plot(x(4,:),'LineWidth',1)
title('x4'); xlabel('Timp')
grid
figure('Position', [550, 550, 900, 600]);
plot(u,'LineWidth',1)
title('Comanda'); xlabel('Timp')
grid

%% Comparatie K variabil si K fix
clc; close all;

figure('Position', [550, 550, 1200, 1200]);
subplot(5,2,1)
plot(time, simulinkOut.simoutResponse.Data(:,1),'LineWidth',1); grid; title('X1 K fix'); xlabel('Timp');
subplot(5,2,2)
plot(x(1,:),'LineWidth',1); grid; title('X1 K variabil'); xlabel('Timp');

subplot(5,2,3)
plot(time, simulinkOut.simoutResponse.Data(:,2),'LineWidth',1); grid; title('X2 K fix'); xlabel('Timp');
subplot(5,2,4)
plot(x(2,:),'LineWidth',1); grid; title('X2 K variabil'); xlabel('Timp');

subplot(5,2,5)
plot(time, simulinkOut.simoutResponse.Data(:,3),'LineWidth',1); grid; title('X3 K fix'); xlabel('Timp');
subplot(5,2,6)
plot(x(3,:),'LineWidth',1); grid; title('X3 K variabil'); xlabel('Timp');

subplot(5,2,7)
plot(time, simulinkOut.simoutResponse.Data(:,4),'LineWidth',1); grid; title('X4 K fix'); xlabel('Timp');
subplot(5,2,8)
plot(x(4,:),'LineWidth',1); grid; title('X4 K variabil'); xlabel('Timp');

subplot(5,2,9)
plot(time, simulinkOut.simoutCommand.Data(:,1),'LineWidth',1); grid; title('Comanda K fix'); xlabel('Timp');
subplot(5,2,10)
plot(u,'LineWidth',1); grid; title('Comanda K variabil'); xlabel('Timp');

%% LQR pentru probleme de urmarire
clc; close all;

Q = eye(4);
R = 1;
C = [0 0 1 0];

K = lqr(A,B,Q,R)
N = -inv(C*inv(A-B*K)*B)
Kz = K(4);
Kx = K(1:3);

simulinkOut = sim('LQR_probleme_de_urmarire.slx');

time = simulinkOut.simoutResponse.Time;
figure('Position', [550, 550, 900, 600]);
subplot(211)
plot(time, simulinkOut.simoutResponse.Data, time, simulinkOut.simoutReference.Data); grid; title('Răspunsul sistemului Q=50*eye(4) R=1');legend('Ieșire', 'Referință');
subplot(212)
plot(time, simulinkOut.simoutCommand.Data); grid; title('Comanda');

figure('Position', [550, 550, 900, 600]);
for i=1:4
    subplot(4,1,i)
    plot(time, simulinkOut.simoutStates.Data(:,i)); grid; 
end

%% LQR pentru probleme de urmarire cu integrator
clc; close all;

C = [0 0 1 0];
Ae = [A zeros(4,1); -C 0];
Be = [B;0];
Q = eye(5);
R = 1;

K = lqr(Ae,Be,Q,R)
Kz = K(5);
Kx = K(1:4);
%% LQR cu estimarea starii
clc; close all;
Q = eye(4);
R = 1;
C = [0 0 1 0];

if rank(obsv(A,C))~=4
    error('Sistemul nu este observabil');
end

K = lqr(A,B,Q,R)
N = -inv(C*inv(A-B*K)*B);

poli = eig( A-B*K)
poli = real(poli)*10;
L = place(A',C',poli)
