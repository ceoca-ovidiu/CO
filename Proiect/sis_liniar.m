function xd=sis_liniar(x,u)
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
     A41 A42 A43 A44]./lambda;

B = [0; B21; 0; B41]./lambda;

xd=A*x+B*u;
end