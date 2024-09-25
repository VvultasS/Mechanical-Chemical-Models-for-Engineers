%% 3-1-8-Fre KBT340.2023.01.10.tenta B4
% clc,clear,
% 9-Tis-19 vardagsrummet% Ska skapa en lokal funktion som ger
% cd,v,Re från kraftbalans(cd,v) och korrelation(cd,Re) 
% och definition(Re,v). 

%Sista frågan på transporttentan hade en förenkling som krävdes
%för att göra problemet algebraiskt lösbart. Men när jag såg den
%lovade jag mig själv att modellera diffekvationen i MATLAB.
%(diametern av droppen minskar under fallet då champangen
%avdunstar, och för att underlätta beräkningen av transport
%ekvationerna gavs ett medelvärde. Jag vill modellera 
%diffekvationen som uppstår. 
%3 okända, 3 samband

% Ons-13:12 Problem now is that there is no acceleration
%in ode so the particle is assumed to be AT terminal velocity,
%which is changing due to evaporation. To Account for this I
%would have to add another element to U. U=[d z v]
format shortg
disp('___________________'),clear; 
% clf
clc
    % Fixed information
Tb=7; Tp=3; p=520; Psat=800; P=101325; 
T=mean([Tp Tb])+273.15; R=8.3145;

    % Momentum Transfer Constants
o.rhoP=999; o.rhoF=1.269; 
o.ny=1.382e-5;
o.vconst=4/3*(o.rhoP/o.rhoF-1);
o.g=@(z) 6.6743e-11*5.972e24/(z+6371e3).^2;

    % Mass Transfer Constants
DABP0=2.634; T1=298.15;  
eAB_k=(97*356)^(1/2);
kT1_e=T1/eAB_k; kT2_e=T/eAB_k;
OmegaD(1)=1.264-(1.264-1.279)/(1.65-1.6)*(1.65-kT1_e);
OmegaD(2)=1.198-(1.198-1.215)/(1.5-1.45)*(1.5-kT2_e);
DAB=DABP0/P*(T/T1)^(3/2)*OmegaD(1)/OmegaD(2);
o.Shconst=.552*(o.ny/DAB)^(1/3);
o.NAconst=-DAB/R/T*log((P-Psat)/(P-p))...
    *18.016e-3/o.rhoP*6;

  
  
    % Input space
z0=324;
v0=1e-5;
D0=[1:-.1:.40];



    % Setup
str1=string(['D0=' num2str(D0(1)) 'mm']);
opt=odeset('Events',@evaporate,'AbsTol',1e-14,...
    'RelTol',1e-11,'Refine',4);
str2=[str1 string(num2cell(D0(2:end)))];
% for j1=1:1,clf(j1),end
Jmax=2000; low=-100;
I=length(D0); J=low*ones(1,I); 
t=low*ones(Jmax,1);
D=low*ones(Jmax,I); Z=D; Dpercent=D; V=D;

    % Model
for i=1:I
    d0=D0(i)*1e-3;
    [t,U,ye,te,ie]=ode15s(@ode,[0 10000],[d0 z0 v0],opt,o);
    J(i)=size(t,1);         T(1:J(i),i)=t;
    d=U(:,1)*1e3;           D(1:J(i),i)=d;
    z=U(:,2);               Z(1:J(i),i)=z;
    v=U(:,3);               V(1:J(i),i)=U(:,3);
    dr=d/d0*100e-3;   Dpercent(1:J(i),i)=dr;
    pt=plot3(v,dr,z); hold on
end

    % Output space
T1=cat(1,T,low*ones(Jmax-max(J),I));
% plot3(T1,Dpercent,z) % I don't remember the purpose of this
xlabel v, ylabel 'd%', zlabel z
title('Champagne Drop'), grid on, box on
legend(str2)
axis([min(v) max(v) 0 100 0 z0-1])
pt.DataTipTemplate.DataTipRows(1).Label = "v";
pt.DataTipTemplate.DataTipRows(2).Label = "d%"; 
pt.DataTipTemplate.DataTipRows(3).Label = "z"; 
dtRows = dataTipTextRow("time",t);
pt.DataTipTemplate.DataTipRows(4) = dtRows;
% dt=datatip(pt,v,dr,z)

    % Functions
function DU=ode(~,U,o), d=U(1); z=U(2); v=U(3);
%     % Momentum
% cd=fzero(@cdsolve,.5,[],U,o);
% v=sqrt(o.vconst*o.g(z)*d/cd);

Re=abs(v*d/o.ny);
    cd=(24./Re)+(2.6.*(Re./5))/(1+(Re./5).^1.52)+...
       (0.411.*(Re./263000).^-7.94)/(1+(Re./263000).^-8);
    a=(1-o.rhoF/o.rhoP)*o.g(z)-3/4/d*v^2*cd/o.rhoP;
    % Mass
Sh=2+o.Shconst*Re^(1/2);
ddot=(o.NAconst*Sh*d)^(1/3);
DU=[-ddot; -v; a]; end

% % obsolete function since ode now tracks acceleration
% function cdOut=cdsolve(cdIn,U,o), d=U(1); z=U(2);
% v= sqrt(o.vconst*o.g(z)*d/cdIn); 
% Re= v.*d./o.ny;
% cdkorr=abs((24./Re)+(2.6.*(Re./5))/(1+(Re./5).^1.52)+...
%        (0.411.*(Re./263000).^-7.94)/(1+(Re./263000).^-8));
% cdOut=cdIn-cdkorr; end

function [R,i,d]=evaporate(~,U,~),R=min(U)-1e-9;i=1;d=0;end
