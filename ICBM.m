%% 3-1-5-Tis-17 
% Trying to do trajectory in polar 
% coordinates
%2023-10-02 3-1-6-Mån-12 Codename: ICBM
% 6-Sön

% 8-Ons The project is scattered and I don't
%feel like implementing more control models for
%the breaking of the fuel system, nor thermodynamical
%functions for the rocket engine. I would need to
%pick a specific fuel (density,specific heat...).
%The rocket geometry is overdefined.
% 3-2-1-Lör No longer the case.

% 3-4-2-Lör-00:52 % Ronnie Andersson recommended that I
%ease into the running-out-of-fuel condition in a milliseconds
%in order to improve ode-solvers convergence. This way it works
%even with ode45, simple and smart.

%%
% 8-Sön-19:57 
clear
disp('__________')
    % External Constants
C.G=6.6743e-11;
C.M=5.972e24;
C.Re=6371e3;
C.Ratm=[0 11 25 87 400 1000]*1e3+C.Re;
C.air= [0       1.225000 101325
        11e3	 .36391  22632.1
        20e3	 .08803  5474.89
        32e3	 .01322  868.02
        47e3	 .001430 110.91
        51e3	 .000860 66.94
        71e3	 .000064 3.96
        100e3   0        0
        100e6   0        0]+[C.Re 0 0];	

    % Rocket
C.ve=3e3; %1.7-4.5km/s | empirical formulae
C.pe=2e5;
cd=.7;
C.Ae=.1;
mfuel=100;
C.m=10;
C.rhoFuel=70.85;
C.rhoP=5000;
C.mdot=1;
rhoMachine=2700; %assume density of aluminium
C.V=C.m/rhoMachine+mfuel/C.rhoFuel;

    % Initial Condition
speed=1e1;
angle=45;


time=[0 6 0; 1 3600 24];
tspan=[0 time(1,:)*flip(cumprod(time(2,:)))'];
opt=odeset('AbsTol',1e-3,'RelTol',1e-3,'Events',@hit);
ag=linspace(-pi,pi,400)';

%         for i=1:2, clf(i),end
% Iterants came from the newer model ChampagneDroppe model

Iterants=[10 30 50 100 150];
% Iterants=[.7 1]; 
I=length(Iterants);
str1=string(['mfuel=' num2str(Iterants(1))]);
% str1=string(['C.Ae=' num2str(Iterants(1))]); % revelation
str2=[str1 string(num2cell(Iterants(2:end)))];
for i=1:I
% C.Ae=Iterants(i);
mfuel=Iterants(i);

C.Ad=cd/2*C.Ae;
vr0=speed*sind(angle);
w0=speed*cosd(angle)/C.Re;
U=[C.Re 0 vr0 w0 mfuel+C.m]';
[t,U,te,ye,ie]=ode45(@Newton2,tspan,U,opt,C);

r=U(:,1); a=U(:,2); vr=U(:,3); w=U(:,4); 
m=U(:,5); wr=w.*r; v=vecnorm([vr wr],2,2);
ratm=r-C.Re; J=size(t,1); ar=a.*r;
    % plots
    figure(1),
    subplot(1,2,1);
    pt=loglog(abs(ar),abs(ratm)); grid on, hold on,
    title('Chartesian Trajectory'); 
    legend(str2);
    top=max(max(ar,ratm));
    axis([0 top 0 top])
    
    pt.DataTipTemplate.DataTipRows(1).Label = "y";
    pt.DataTipTemplate.DataTipRows(2).Label = "x"; 
    dtRows = [dataTipTextRow("time",t),dataTipTextRow("mass",m)];
    pt.DataTipTemplate.DataTipRows(3:4) = dtRows;
% if max(abs(a))>.1    
    axis([-.1 ar(end)+.1 0 max(ratm)+1])
    subplot(1,2,2);
    polarplot(a,r,'LineWidth',1), hold on
    title('Trajectory')
% end
% % polarplot(ag,C.Ratm.*ones(size(ag)),...
% %     'LineStyle','--','color','k'),hold on,
% %  legend('Exo','Thermo','Meso','Strato','Tropo',...
% %      'Surface','Location','best')
% 
% polarplot(ag,C.Ratm(1).*ones(size(ag)),...
%     'LineStyle','--','color','k'),
% polarplot(ag,C.Ratm(end).*ones(size(ag))...
%     ,'LineStyle','--','color','b')
% legend('Surface','Exo','1','2','Location','best')
% if J<100, for j=1:J-1
%         polarscatter(a(j),r(j),'k','.');
%     end ,else, polarplot(a,r,'k'), end,
%  legend('Surface','Tropo','Strato','Meso',...
%      'Thermo','Exo','Location','best')
% 

figure(2), 
 subplot(2,2,1),title('Kinematics')
  loglog((t),(v)), hold on, grid on,
%   loglog((t),abs(wr)),
%   legend('speed','tangential','Location','best'),
  ylabel('speed'), 
  
 subplot(2,2,2)
  loglog(ratm,v), hold on, grid on
%   loglog((ratm),abs(wr)),
%   legend('speed','tangential','Location','best'),
  
  subplot(2,2,3),
  plot(log10(t),(m-C.m)/mfuel), hold on, grid on,
  xlabel('log10 seconds'), ylabel('% fuel'),

 subplot(2,2,4)
  plot(log10(abs(ratm)),abs((m-C.m)/mfuel)), 
  hold on, grid on,
  xlabel('log10 meters above surface'), 


%Hör till fig 1
% for j1=1:6, yline(C.Ratm(j1)-C.Re),end
%  legend('Exo','Thermo','Meso','Strato','Tropo',...
%      'Surface','Location','best')

 apex=max(r)-C.Re;
end
% -----------------------------------
function DU=Newton2(~,U,C)
r=U(1); vr=U(3); w=U(4); m=U(5);
wr=w*r;  v=norm([vr wr]); 
i=sum(r>=C.air(:,1)); %Elapsed time is 0.038358 seconds.
if i==0, i=1;end
rhoF=C.air(i,2);
p0=C.air(i,3);

rhoB=m/C.V;

% THRUSTER CONTROLLER
if m-C.m>0
mdot=C.mdot;
T=mdot.*C.ve + (C.pe-p0).*C.Ae;
% T=mdot*C.ve;
elseif m-.9*C.m>0
mdot=C.mdot*(1-C.m/m);
T=mdot.*C.ve;
else, mdot=0; T=0;
end


boy=(1-rhoF/rhoB);
    %acceleration components
g=C.G*C.M*boy/r^2;
c=w*wr;
D=C.Ad*rhoF;

DU=[vr w [c-g 0]+[vr w]*(T/v-D*v)/m -mdot]';
end
function [R,i,d]=hit(~,U,C)
R=min((U(1)-C.Re+.5),1.4*C.Ratm(end)-U(1)); %hitting the ground or blasting off to the far beyond
i=1; d=0;
end


% clear all
% % C.Ratm=[0 11 25 87 400 1000]*1e3+C.Re;
% C.air= [  0   1.225000 101325
%         11e3	 .36391  22632.1
%         20e3	 .08803  5474.89
%         32e3	 .01322  868.02
%         47e3	 .001430 110.91
%         51e3	 .000860 66.94
%         71e3	 .000064 3.96
%         100e3 0        0
% 				];
% x=C.air(:,1); y=C.air(:,2);		
% clc
% r=0:1:51e3;
% J=length(r);
% rhoF=zeros(J,1);
% tic
% for j=1:J
% rhoF(j)=interp1(x,y,r(j),'previous'); Elapsed time is 0.316305 seconds.
% % i=sum(r(j)>=C.air(:,1)); %Elapsed time is 0.038358 seconds.
% % rhoF1(j)=C.air(i,2);
% end
% toc
