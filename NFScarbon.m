%% 2023-03-10  2-3-8-Fre
% NFS Carbon

% when the car is at coordinate along the road x in the 
% intervall (x1,x2) its maximal velocity is vlimit12, meaning it must
% deaccelerate to vlimit12. The code is purely kinematic to begin with.
% Idea: Define its a a and calculate its velocity v within x12
% by integrating v= a dt.  and obtain t  through  v = dx/dt => 
% dt = 1/v dx;  
% Take an input of an inital velocity v0 and an a vector A. 
% Render an output of time t
% v = dx/dt => dx = v dt
% a = acceleration, v= velocity , x= coordinate along road
% a = dv/dt =>  a = dv/dt * v/(dx/dt) = v dv/dx
% dt= 1/v dx
% a dx = v^2/2 => v= sqrt(2 a dx)
%% input V X n % First generation
% % % % clear
% % % % clc
% % % % v0=1;
% % % % L=[35 10 40 5 100];
% % % % a=[10 5 30 10 100];
% % % % N=100; M=5;
% % % % t=zeros(N,M); v=t;
% % % % v(1)=v0;
% % % % dx=L/N;
% % % % for m=1:M
% % % %     for n=1:N
% % % %         v(n+1,m) = sqrt(2*dx(m)*a(m)) + v(n,m);
% % % %         t(n+1,m) = dx(m)./v(n,m);
% % % %     end
% % % %     if m==M
% % % %         break
% % % %     end
% % % %     v(1,m+1)=v(end,m);
% % % % end
% % % % v
% % % % timeToClearEachCorner=sum(t,1)
% % % % 
% % % % Ttot=sum(t,'all')
%% WORKING acceleration model % Second generation

% acceleration is going to be a function of vlimit a couple
% dx steps ahead. 
% a = exp(vlimit-v)

% close
clear
clc

M=30;
random1=rand(1,M);
random2=rand(1,M);

mass=1200; Pmax=200e3;  fStat=.8; 
N=100;
drag=1.204*.30*.7/2;
gravity=9.81; v0=20;

%
L= 200.*random1;


Vlimit=sqrt(fStat*gravity*random2*200);
T=zeros(N,M); V=T; A=T; X=T;
V(1)=v0;


Dx=L/N;                         % Too be replaced

thrAggro=6; brkAggro=11;
crash=.5;
for m=1:M
    for n=1:N
        v=V(n,m);
        dx=Dx(m);
        vlimit=Vlimit(m);
%         if A(n,m)< s
%             friction =
%         end
%         vlimit=sqrt(friction*gravity*random2*1000);
        dv=(v-vlimit)./vlimit;
        if     dv<0
            thr= 2.*abs(1./(1+exp(dv*thrAggro))-.5);
            brk=0;
        elseif dv>0 
            thr=0;
            brk= 2.*abs(1./(1+exp(dv*brkAggro))-.5);
        elseif dv>crash
            disp('crash')
            break
        else 
            thr=0;
            brk=0;
        end
        A(n+1,m) =(thr*Pmax/v - fStat*mass*gravity*brk -(drag*v^2 + 0.005+(0.01+0.0095.*(v./100.*3.6).^2)./2.9))/mass;
        V(n+1,m) = sqrt(2*dx*A(n+1,m) + v^2);
        T(n+1,m) = dx/v + T(n,m);
        Thr(n,m)=thr;
        Brk(n,m)=brk;
        X(n+1,m)=dx+X(n,m);
    end
    if m==M
        break
    end
    V(1,m+1)=V(end,m);
    T(1,m+1)=T(end,m);
    X(1,m+1)=X(end,m);
end
vlin=V(:);
tlin=T(:);
xlin=X(:);
Ttotal=T(end,end);

%%
figure(1)
plot(tlin,xlin)
%%
% close all
figure(3)
xp=X(end,:);
vp=Vlimit;
repxp=repmat(xp,2,1); rx=repxp(:)'; rx=[0 rx(1:end-1)];
repvp=repmat(vp,2,1); rv=repvp(:)';
hold on
plot(xlin,vlin,'black','LineWidth',2)
plot(rx,rv,'blue','LineWidth',1)

% Make a function a(v) such that if v large then a small,
% and if

% a(F)= F/m; F(P)=P/v - Resistance(v); 
% P(Throttle)=Throttle*Pmax; Throttle=Real(0,1);

% 00:20 I really need to sleep now.
%
%% throttle = fcn that takes v and vlimit and returns when v-vlimit small 
% then throttle is small, and if v-vlimit
% also larger v values will give larger throttle, as the engine struggles
% against the added friction brought on by greater v.
% if vlimit=infinite => 
% plot(v,throttle)

%% 00:55 Goodnight

%% Demonstration of the sigmoid-esque control function
thrAggro=4;
brkAggro=1000;
dvt=-1:.01:0;


thr= 2.*abs(1./(1+exp(dvt*thrAggro))-.5);

dvb=0:.01:1;
brk= 2.*abs(1./(1+exp(dvb*brkAggro))-.5);
dv=[dvt,dvb];
control=[thr brk];
figure(4)
hold on 
plot(dv,abs(dv))
plot(dv,control)
% yline(0.5)
% xline(0)
hold off

%% Vlimit
% manufacture Vlimit as a function of something more road-specific
% friction*gravity == vlimit^2/r
clc
friction=.5;
coeffs=1e3.*rand(10,1);
r=abs(sum(coeffs.*sin(L)));
Vlimit=sqrt(friction*gravity*r);

%% Varying dx 
% plot(xp,vlimit)
X=[];
dx=[1 2 3];
% x=[0:dx(1):L(1) L(1):dx(2):L(2)]';
for i=1:3
   x=L(i):dx(i):L(i+1);
   X=[X x];
end

%%
figure
f=x.^2-1000.*cos(x);
scatter(x,f)
diff(f)

