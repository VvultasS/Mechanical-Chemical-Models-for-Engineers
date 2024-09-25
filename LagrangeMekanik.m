%% 3-3-4-Fre-20:45 Lagrangemekanik     
%kolla semilogx!

clear
clc
close all
format short g
tic
global g l km Im Dm
m=1; g=-9.81; k=800; l=1; 
I=1 ; %internt motstånd
% D=10; %drag 
tspan=15;
initial=[0 90 0 0]+[0 180 0 0];
o=odeset('AbsTol',1e-8,'RelTol',1e-8,'Refine',4);
s=[0 1 10]; %stryrvariabel
% sStorhet='k= '; sEnhet=' N/m';
sStorhet='D= '; sEnhet=[];
str1=string([sStorhet num2str(s(1))]);
str2=string([num2str(s(end)) sEnhet]);
str3=[str1 string(num2cell(s(2:end-1))) str2];
legendstring=[str3 "initial" "Jämviktsradie"];
ystring=["x" "h" "xd" "hd"];
for j=1:numel(s)
  D=s(j);
km=k/m; Im=I/m; Dm=D/m*pi/180;
[t,u]=ode15s(@L,[0 tspan],initial,o); 
t=t.';u=u.';
x=u(1,:); h=u(2,:); xd=u(3,:); hd=u(4,:);
toc

fig2=figure(2);
for i=1:4
	subplot(2,2,i),
	plot(t,u(i,:)), hold on
	xlabel t
	ylabel(ystring(i))
  if i==1
    [px,p]=max(abs(u(i,:))); 
% p=floor(numel(t)*.7);
% px=u(i,p);
    text(t(p),px+6e-8,char(str3(j))),
  end
end

fig3=figure(3);
r=x+l; J=numel(t)
hl=linspace(0,360,J);
xr=r.*sind(h);
yr=r.*cosd(h); 
plot(xr,yr),hold on
% polarplot( r+l,h*pi/180),hold on
% polarplot(r,h)


end
scatter(xr(1),yr(1),'filled','k')
xl=l.*sind(hl);
yl=l.*cosd(hl);
plot(xl,yl)
legend(legendstring)
fig2.Position=([62 558 1056 798]);
fig3.Position=([201 -3 689 501]);
figure(3),legend(legendstring)

%%
figure(3)
legend off
for j=1:80:J-1
	scatter(xr(j),yr(j),'.','k')
	hold on
	drawnow limitrate
end

function du=L(~,u), global g l km Im Dm
x=u(1);h=u(2);xd=u(3);hd=u(4);
xdd=+(l+x)*hd^2    +g*cosd(h) -km*x -Im*xd;
hdd=-2*xd*hd/(l+x) -g*sind(h) -Dm*(x+l)*hd;
du=[xd;hd;xdd;hdd]; end
