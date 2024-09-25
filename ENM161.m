%% 3-2-8-Tis-11 till 22
clear
% clc
disp _____________
m4=192.5;
format bank
   %H2 CO  CO2 CH4
xb=[.4 .25 .25 .1]';
    %N2     O2    CO2   H2O
wl=[.7573 .1743 .0376 .0308]';
Ml=[28.04      32     44.02   18.02]';
B=136.68e-3;
XL=wl./Ml;
XLsum=sum(XL)
xl=XL/XLsum; xlans=xl*100
% H2 + .5O2 -> H2O
v=[  0   0  0  0 %N2
   -.5 -.5  0 -2 %O2
     0   1  1  1 %CO2
     1   0  0  2];%H2O
L=m4*XLsum;
v*xb
format short
R= xl*L + v*xb*B
format bank
xr=R/sum(R)*100
%%
p=1;
n=67;
k=(n-1)/2/p





%%
clc
clear
format shortg

L=100;
B=10;
o.xL=[.79 .21 .00 .00 .00]';
o.FL=L*o.xL;
xb=[.1 .9]';
o.FB=xb*B;
v=[0 0; -1.5 -1.25; 1 1; 2 1.5; 0 1];
o.FR=o.FL + v*o.FB;
R=sum(o.FR);
xR=o.FR/R
% EnergyBalance
% cp should be a fcn of T1,T2. From integral of a
% polynomial with coeffs living in rank 2; accepting
% scalar input, and returning a vector that lives in 
% rank 1.
%funktioner
cap=[28.883  -.157  .808 -2.871  %N2
     25.460  1.519 -.715   1.311 %O2
     22.2243 5.977 -3.499  7.464 %CO2
     32.218  .192  1.055 -3.593  %H2O
     29.35   -.094  .9747 -4.187 %NO
].* [1       1e-2   1e-5   1e-9];
n=[1 2 3 4]';
cp=@(T1,T2)cap*((T2.^n - T1.^n)./n);

%input
o.D=2;
o.L=1;
o.T1=300;
o.Tr=298;
o.P1=1;
o.P2=15;
o.etaKis=.9;
o.etaTis=.85;
o.R=8.3145;
T2=667.5;
o.M=[28 32 44 18 30]';
for i=1:4
    cp12=o.xL'*cp(o.T1,T2)./(T2-o.T1);
    kappa=cp12/(cp12-o.R);
    tau=(kappa-1)/kappa;
    T2is=o.T1/(o.P1/o.P2)^tau;
    T2=o.T1+(T2is-o.T1)/o.etaKis;
end
o.T2=T2;

%beräkning
o.A=pi*o.D*o.L;
%alla l är tagna från WRF vid Tf=1000K
l.beta=.803e6*o.D^3; 
l.k=0;% l.k=6.7544e-2;
l.Pr=.702;
o.H=[725.7e3 82.05e3]';
Q.Q2=o.FL'*cp(o.Tr,T2);
Q.Qbr=o.FB'*o.H;
T3=fzero(@(T3)BrKammare(T3,l,o,Q,cp),1106.8)
[~,Q]=BrKammare(T3,l,o,Q,cp);
disp(Q)
Tm=(o.T2+T3)/2;
Tf=(Tm+o.T1)/2;

T4=Tf;
for i=1:4
    cp34=xR'*cp(T3,T4)./(T4-T3);
    kappa=cp34/(cp34-o.R);
    tau=(kappa-1)/kappa;
    T4is=T3/(o.P2/o.P1)^tau;
    T4=T3+(T4is-T3)*o.etaTis;
end
etaMG=.95;
Wk=L*cp12*(T2-o.T1)/etaMG;
WT=R*cp34*(T3-T4)*etaMG;
Wel=(WT-Wk);
Q4=o.FR'*cp(o.Tr,T4)
etaTOT=(Wel+Q4)/Q.Qbr
alfa=Wel/Q4

function [EB,Q]=BrKammare(T3,l,o,Q,cp)
Tm=(o.T2+T3)/2;
% Tf=(Tm+o.T1)/2; 
dT=Tm-o.T1;
Gr=l.beta*dT;
Ra=Gr*l.Pr;
if Ra<1e-5
    error('Ra<1e-5')
elseif Ra>1e12
    error('Ra>1e12')
end
Nu=(.6+.387*Ra^(1/6)/(1+(.559/l.Pr)^(9/16))^(8/27))^2;
h=l.k/o.D*Nu;
Q.Qh=h*o.A*dT;

Q.Q3=o.FR'*cp(o.Tr,T3);
EB=Q.Q2 + Q.Qbr - Q.Q3 - Q.Qh;
end
