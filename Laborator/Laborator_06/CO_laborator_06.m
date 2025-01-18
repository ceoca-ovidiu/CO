clc; close all; clear all;
dt=0.001; %;pas de integrare
t=0:dt:2;
N=length(t); %nr puncte

x=1;
u=ones(1,N);
eps=1e-4;
lambda(N)=0;
dH=1;
i=1;
s=0.01;

while norm(dH)>eps && i<1e5
    for k=1:N-1
        K1=ec_stare(x(:,k),u(k));
        K2=ec_stare(x(:,k)+0.5*dt*K1,u(k));
        K3=ec_stare(x(:,k)+0.5*dt*K2,u(k));
        K4=ec_stare(x(:,k)+dt*K3,u(k));
        
        x(:,k+1)=x(:,k)+1/6*dt*(K1+2*K2+2*K3+K4);
    end
    
    for k=N:-1:2
        K1=ec_lambda(lambda(k));
        K2=ec_lambda(lambda(k)-0.5*dt*K1);
        K3=ec_lambda(lambda(k)-0.5*dt*K2);
        K4=ec_lambda(lambda(k)-dt*K3);
        
        lambda(k-1)=lambda(k)-1/6*dt*(K1+2*K2+2*K3+K4);
    end
    i=i+1;
    dH=2*u+lambda;
    if i>13000
        s=1e-4;
    end
    u=u-s*(2*u+lambda)/norm(2*u+lambda);
     norm(2*u+lambda)
end
figure;
plot(x); grid;
figure;
plot(u); grid;