%%2024aug27
clc
clear
global S ub
Tb=30;
yb=0;
ym=.74;
ub=cat(1,Tb,yb);
L=3e-2;
d=5e-3;
R=13e-2;
PWR=3;
S_A=PWR/(d*pi*R^2);
S=@(t)S_A*cos(t/5)^2;
N1=10;
N2=200;
N=N1+N2;
h1=d/N1;
x1=h1:h1:d;
h2=(L-d)/N2;
x2=d+h2:h2:L;
x=[x1 x2];
steak=logical([zeros(1,N1) ones(1,N2)]);
pan=logical(1-steak);

M1=50;
M2=M1;
M3=M1;
t0=linspace(0,100,M1);
t1=linspace(2*t0(end)-t0(end-1),200,M2);
t2=linspace(2*t1(end)-t1(end-1),300,M3);
t=[t0 t1 t2];

u0i=[Tb;0];
IC0=@(x) u0i;
u0pan=pdepe(0,@PDE,IC0,@BC,x1,t0);

u0=[u0pan cat(3,Tb,ym).*ones(M1,N2)];
u1i=permute(u0(end,:,:),[3 2 1]);
IC1=@(X) u1i(:,X==x);

u1=pdepe(0,@PDE,IC1,@BC,x,t1);

u2i=permute(u1(end,:,:),[3 2 1]);

u2if(:,steak)=flip(u2i(:,steak),2);
u2if(:,pan)=u2i(:,pan);

IC2=@(X) u2if(:,X==x);

u2=pdepe(0,@PDE,IC2,@BC,x,t2);

clf

T0pan=u0pan(:,:,1);

T1=u1(:,:,1);
T2=u2(:,:,1);
y1=u1(:,:,2);
y2=u2(:,:,2);
%
figure(1)
clf
mesh(t0,x1,T0pan')
hold on
mesh(t1,x,T1')
mesh(t2,x,T2')
title 'Temperature Distribution'
xlabel seconds
ylabel meters
zlabel celcius

figure(2)
clf
mesh(t1,x,y1')
hold on
mesh(t2,x,y2')
title 'Water Mass Fraction Distribution'
xlabel 'seconds'
ylabel 'meters'

function [c,f,s]=PDE(x,t,u,dudx)
global S 

if x<.0050
  c=[3501 ; 1];
  f=[71 0 ; 0 0]*dudx;
  s=[S(t) ; 0];
elseif x>.0055
  c=[1e4 ; 1];
  f=[1.23e-3 0 ; 0 1e-8*u(1)]*dudx;
  s=[0 0 ; 0 0]*u;
else
  c=[3501 ; 1];
  f=[1e-2 0 ; 0 1e-7*u(1)]*dudx;
  s=[-20 0 ; 0 .01]*(u.^2);
end
end

function [pl,ql,pr,qr]=BC(xl,ul,xr,ur,~)
% p + q*f == 0
global ub

pl=[0;0];
ql=[1;1];
% pl=cat(3,0,0);
% ql=cat(3,1,1);

h=.5;
kc=1e-3;
pr=[h 0;0 kc]*(ur-ub);
qr=[1;1];
end



% if-sats        -punktvis i PDE 
% konkatenation  -punktvis i IC 
% event 


% antingen testa event igen
% eller transformera vektorn u0
% iflip=logical([zeros(1,N1) ones(1,N2)]);
% u0(logical(1-iflip))=u20(logical(1-iflip))
% u0(iflip)=flip(u20(iflip)');
% till en skalärfunktion av x.
