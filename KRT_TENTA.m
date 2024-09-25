clear
clc
%På teoriuppgiften 4.f) tentan 2023-05-29 i kemisk
%reaktionsteknik, skulle man skriva ner ekvationerna 
%som behövdes för att lösa en andra ordningens
%vätskefasrxn i tubreaktor med filmtransport.

%För att säkerställa poäng skrev jag en
%MATLAB-rutin som kunde lösa problemet.
%Kvällen därpå replikerade jag, från minne, denna rutin.
%Den fungerade som förväntat.
Sh=2;
dp=1e-3;
Da=1e-5;
sp=6/dp;
kc=Sh*Da/dp*sp;
k=0.05;
opt1.AbsTol=1e-12; opt1.Refine=4;
L=[0 10];
caf=400;
 t=zeros(10,2);
for j=1:10
tic
[~,cab1]=ode45(@kinetik1,L,caf,opt1,kc,k);
t(j,1)=toc;
%uppdaterad snabbare, enklare, precisare algoritm
tic
[z,cab2]=ode45(@kinetik2,L,caf,opt1,kc,k);
t(j,2)=toc;
end
meantime=mean(t,1)
meandiff=mean(cab1-cab2)
figure
plot(z,cab1), grid on, hold on, axis([L 0 caf])

%%
function dcabdz=kinetik1(~,cab,kc,k)
    cas=fsolve(@find_cas,cab,optimset('display','off'));
    dcabdz=-2*k*cas^2;
    function res=find_cas(cas)
        res=2*k*cas^2-kc*(cab-cas);
    end   
end
function dcabdz=kinetik2(~,cab,kc,k) 
    cas=fzero(@(cas)2*k*cas^2-kc*(cab-cas),cab);
    dcabdz=-2*k*cas^2;
end
