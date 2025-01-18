clc; close all; clear all;

syms x(t) l(t)

Dx=diff(x);
Dl=diff(l);

s=dsolve(Dx==x-l/2,Dl==-2*x-l,x(0)==1,l(5)==0)

subplot(121)
fplot(s.l); grid;
subplot(122)
fplot(s.x); grid;

%% masinuta
clc; close all; clear all;

syms x(t) l(t)

Dx=diff(x);
Dl=diff(l);
a1 = 1;
a2 = 1;

s=dsolve(Dx==-2*x-l/(2*a1),Dl==2*a2*(x-cos(t))-2*l,x(0)==0,l(20)==0)

subplot(121)
fplot(s.l); grid; 
subplot(122)
fplot(s.x); grid; 