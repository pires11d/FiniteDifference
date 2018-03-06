clc
clf
clear all
close all

n = 10;
mf = n;
nf = n;
L = 5;
R = 1;
dz = L/(n-1);
dr = R/(n-1);
dt = 0.05;
tf = 5;

Dz = 0.001;
Dr = 0.01;
k = 1;
Ka = 0;
u = 1.5;

z = 0:dz:L;
r = 0:dr:R;
t = 0:dt:tf;

pf = length(t);

Ca = zeros(mf,nf,pf);
Cao = 0;
Cae = 100;

Camin = min(Cao,Cae);
Camax = max(Cao,Cae);
Lmax = max(R,L);

Ca(:,:,1) = Cao;        %Ca(z,r,0) 
Ca(1,:,:) = Cae;        %Ca(0,r,t)

for p=1:pf-1
    for m = 2:mf-1
        for n = 2:nf-1
            Ca(m,1,p)=Ca(m,2,p);
            adveccao_z = -u*(Ca(m,n,p)-Ca(m-1,n,p))/(dz);
            difusao_z = Dz*(Ca(m+1,n,p)-2*Ca(m,n,p)+Ca(m-1,n,p))/(dz^2);
            difusao_r = Dr*((1/((n-1)*(dr^2)))*(Ca(m,n+1,p)-Ca(m,n,p))+(Ca(m,n+1,p)-2*Ca(m,n,p)+Ca(m,n-1,p))/(dr^2));
            reacao = -k*Ca(m,n,p)/(1+Ka*Ca(m,n,p));
            Ca(m,n,p+1) = Ca(m,n,p) + dt*(adveccao_z + difusao_z + difusao_r + reacao);
        end
    end
    Ca(mf,:,p)=2*Ca(mf-1,:,p)-Ca(mf-2,:,p);
    Ca(:,nf,p)=Ca(:,nf-1,p);
end

for p=1:pf-1
figure(1)
surf(r,z,Ca(:,:,p))
hold on
title('Concentração de A no Reator','FontSize',14)
surf(-r,z,Ca(:,:,p))
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Ca (mol/L)')
caxis([Camin,Camax])
view(-140,20)
hold off
xlim([-Lmax Lmax])
ylim([0 Lmax])
zlim([Camin Camax])
drawnow
end



