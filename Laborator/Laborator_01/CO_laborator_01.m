clc; close all; clear all;

syms x1(t) x2(t)
 
s = dsolve(diff(x1, t) == x2, diff(x2, t) == -x1, x1(0)==0, x2(0)==1)

figure;
fplot(s.x1, [0, 3*pi]); hold on;
fplot(s.x2, [0, 3*pi]);

%% facem cu ode acum 
close all; clc;

[t,x] = ode23(@xprime, [0, 3*pi], [0 ; 1]); % [0;1] conditii initiale
                                            % pot pune in loc de [0, 3*pi]
                                            % si 0:0.1:3*pi

figure;
plot(t, x(:,1), t, x(:,2)); 
 
%% Metoda Euler
clc; close all; clear all;

h = 0.01; % cu 0.1 merge spre instabilitate, dar cu 0.01 e mai ok 
tf = 3*pi;
t = 0:h:tf;
N = length(t);
x = [0;1];
final_value = [;];

for i=1:N
    x = x + h*xprime(0.1,x);
    final_value = cat(2,final_value,x);
end

figure;
plot(t,final_value(1,:),t,final_value(2,:))

%% Metoda Runge-Kutta de ordin 4
clc; close all; clear all;

h = 0.1; 
tf = 3 * pi;
t = 0:h:tf;
N = length(t);
x = [0 ; 1];
final_value = [;];

for i=1:N
    K1 = xprime(0.1, x);
    K2 = xprime(0.1, x+(h*K1)/2);
    K3 = xprime(0.1, x+(h*K2)/2);
    K4 = xprime(0.1, x+h*K3);

    x = x + (h*(K1 + 2*K2 + 2*K3 + K4))/6 ;

    final_value = cat(2, final_value, x);
end

figure;
plot(t, final_value(1,:), t, final_value(2,:))

%% Integram de la capat la cap daca am conditii finale nu initiale (Euler intors)
clc; close all; clear all;

h = 0.01;
tf = 3*pi;
t = 0:h:tf;
N = length(t);
x = [-1;0];
final_value = [;];

for i=1:N
    x = x - h*xprime(0.1,x);
    final_value = cat(2,x,final_value);
end

figure;
plot(t,final_value(1,:),t,final_value(2,:))

%% Exercitiu molocotiva (Simulink)
M1 = 1; 
M2 = 0.5;
k = 1;
miu = 0.002;
g = 9.8;

z1 = x1
z2 = diff(x1)

