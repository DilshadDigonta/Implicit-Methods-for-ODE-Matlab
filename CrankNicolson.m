function CrankNicolson()
% Program MP3p02.m
 
clear all

s = pi/18;
% Solving the ODE for
% 0<t<4 with x1(0)=s and x2(0)=0

%Own Crank Nicolson Solver

[t1,x1]=MyCrankNicolson(@dxdt,[0 4],[s 0],0.2);
[t2,x2]=MyCrankNicolson(@dxdt,[0 4],[s 0],0.1);
[t3,x3]=MyCrankNicolson(@dxdt,[0 4],[s 0],0.05);
%Matlab Solver
[tmat,xmat]=ode23(@dxdt,[0 4],[s 0]);
plot(t1,x1(:,1),'k-',t2,x2(:,1),'b-',t3,x3(:,1),'g-',tmat,xmat(:,1),'ro-')
legend('h=0.2','h=0.1','h=0.05','ode23()');
xlabel('t');
ylabel('x');
title("CrankNicolson");



function [t,x]=MyCrankNicolson(ODEfunc,tspan,x0,h)
t=tspan(1):h:tspan(2);
x(1,:)=x0;
for n=1:length(t)-1
xn=x(n,:);
tn=t(n);
tnp1=t(n+1);
x(n+1,:)=fsolve(@(xnp1)CrankNicolsonODEfunc(ODEfunc,tn,tnp1,xn,xnp1,h),x(n,:));
end
function residual=CrankNicolsonODEfunc(ODE,tn,tnp1,xn,xnp1,h)
residual=(xn+(h/2)*(ODE(tnp1,xnp1)+ODE(tn,xn))'-xnp1);


function xp=dxdt(t,x)
xp(1)=x(2);
xp(2)=-16.35.*x(1);
xp=xp';
