clc; close all; clear all;
u_min = -1;
u_max = 1;
x_min = 0;
x_max = 3/2;
pace = 0.5;

u = u_min:pace:u_max;
x = x_min:pace:x_max;

JN = [];
N = 3;

J_star = getH(x);
J_aux = [];
J = [];
MUoptim =[];


for k=N-1:-1:1
    for i=1:length(x)
        for j=1:length(u)
            x_urm = getSum(x(i),u(j));
            if x_urm >= x_min & x_urm <= x_max
                J_aux(j) = interp1(x, J_star, x_urm) + getL(u(j));
            else
                J_aux(j) = 999999;
            end
        end
        [V, I] = min(J_aux);
        J(i) = V(1);
        MUoptim(i, k) = u(I(1));
    end
    J_star = J;
end

MUoptim

xo = [];
uoptim = [];

xo(1) = 1;

for k=1:N-1
    uoptim(k) = interp1(x, MUoptim(:,k), xo(k));
    xo(k+1) = getSum(xo(k), uoptim(k));
end

plot(xo); grid; hold on;
plot(uoptim); grid;

%%
clc; close all; clear all;

Hf = tf(10,[1 10])








