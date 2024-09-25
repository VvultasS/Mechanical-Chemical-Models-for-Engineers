%% 3-4-3-LÖR-17:15 2h
clear
disp __________
global a L R
a=45;
L=1;
R=1e0;
t0=[0 10];
v0=1e-2;
s0=0;
U0=[0 v0 s0];
o=odeset('Events',@fall);
[t,U]=ode15s(@TEMPUS,t0,U0,o);
x=U(:,1); v=U(:,2); s=U(:,3);

time=t(end)
top_speed=max(v)

figure(5)
subplot(3,1,1)
plot(t,x); title x, ylabel [m], hold on
subplot(3,1,2)
plot(t,v); title v, ylabel [m/s], hold on
subplot(3,1,3)
plot(t,s); title slip, ylabel [ ], hold on
xlabel t
function Udot=TEMPUS(~,U)
global a R
v=U(2); s=U(3);
sc=(1-s);
FKMAX=10;
FK=FKMAX*(1-s^2);
% if s<1e-10,  FK=FKMAX; else, FK=FKMAX*(1-s^2); end
% FK=1000;

RR=0;
RMI=1/R;
g=9.82;
vdot=g*(cosd(a)-sind(a)*(s*FK+sc*RR));
sdot=sc*vdot/v - RMI*s*FK*g*sind(a);
Udot=[v;vdot;sdot];
end

function [R,I,d]=fall(~,U)
global L
R=L-U(1);
I=1;
d=0;
end
