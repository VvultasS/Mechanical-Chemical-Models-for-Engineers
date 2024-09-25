 %% 2023-05-11 2-4-8-Tors

clc
clear

%kula input
rhoB=7000;
D=1e-1;
initialHeight=1;
initialSpeed=100;
angle=15;
tfinal=17;


%omgivning
WaterLevel=0;
g=[0 0 -9.82];
% vind=randn(10,3);
vind=[0 5 0]';
%Beräkning
v0=initialSpeed*[cosd(angle) 0 sind(angle)];
r0=[0 0 initialHeight];
geo= 3/4/D/rhoB; % geometry coeff Ap/m/2 
% Ap=pi/4*D^2; % m=pi/6*D^3*rhoB;

% opt=odeset('Events',@hit);
% opt=odeset('refine',4,'AbsTol',1e-1);
opt=[];

[t,U]=ode45(@AirResistance,[0 tfinal],[v0,r0]'...
    ,opt,geo,D,g,vind,WaterLevel);

r=U(:,4:6); x=r(:,1);y=r(:,2);z=r(:,3);
figure(4),clf
plot3(x,y,z,'Color','k','LineWidth',2); hold on
[xMesh,yMesh]=meshgrid(x,y);
watersurf=surf(xMesh,yMesh,zeros(size(xMesh)),...
'EdgeColor','none','FaceColor','b','FaceAlpha',.4);
hold on, box on, grid on
[~,ind]=min(abs(z));
% figure(5)
% axis([0 x(ind) -1 max(y)+1 0 max(z)])
% Re = logspace(-1, 6, 1000); 
% loglog(Re, cd(Re));  % Plotting the correlation in a log-log scale
% xlabel('Reynolds Number (Re)');
% ylabel('Drag Coefficient (Cd)');
% title('Drag Coefficient of a Sphere vs. Reynolds Number');
function dUdt=AirResistance(~,U,geo,D,g,vind,WaterLevel)
% v=U(1:3)'-vind(randi(height(vind)),:)';
v=U(1:3)'-vind';
if U(6)<WaterLevel, ny=.995e-6; rhoF=999;
else, ny=1.56e-5;  rhoF=1.204; 
end
Re= D.*norm(v)./ny;
cd = 24./Re .* (1 + 0.15.*Re.^0.687) ...
    +(2.6.*Re./5)./(1+(Re./5).^1.52)+...
    (0.411.*(Re./263000).^-7.94)./(1+...
    (Re./263000).^-8);
ad=rhoF.*v.^2.*cd.*geo;
a=-sign(v).*ad + g; 
dUdt=[a v]';
end
function [R,isterminal,direction]=hit(t,U)
R=U(6)-100;
isterminal=1;
direction=0;
end
