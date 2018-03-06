clc
clf
clear all
close all

n = 10;
mf = n;
nf = n;
L = 1.0;
R = 0.4;
dz = L/(n-1);
dr = R/(n-1);
dt = 0.1;
tf = 20;

D = 0.005;
k = 1;

z = 0:dz:L;
r = 0:dr:R;
t = 0:dt:tf;

pf = length(t);

Ca = zeros(mf,nf,pf);
Cao = 100;
Cab = 0;

Camin = min(Cao,Cab);
Camax = max(Cao,Cab);
Lmax = max(R,L);

Ca(:,:,1) = Cao;        %Ca(z,r,0) 
Ca(mf,:,:) = Cab;       %Ca(0,r,t)
Ca(:,nf,:) = Cab;       %Ca(z,R,t)

for p=1:pf-1
    for m = 2:mf-1
        for n = 2:nf-1
            Ca(m,1,p+1)=Ca(m,2,p+1);
            Ca(1,n,p+1)=Ca(2,n,p+1);
            Ca(m,n,p+1) = Ca(m,n,p) + dt*D*((Ca(m+1,n,p)-2*Ca(m,n,p)+Ca(m-1,n,p))/(dz^2)+(2/((n-1)*(dr^2)))*(Ca(m,n+1,p)-Ca(m,n,p))+(Ca(m,n+1,p)-2*Ca(m,n,p)+Ca(m,n-1,p))/(dr^2));
        end
    end
    Ca(1,1,p+1)=(Ca(2,1,p+1)+Ca(1,2,p+1))/2;
end

for p=1:pf-1
figure(1)
surf(r,z,Ca(:,:,p))
hold on
surf(-r,z,Ca(:,:,p))
surf(r,-z,Ca(:,:,p))
surf(-r,-z,Ca(:,:,p))
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Concentracao (mol/L)')
caxis([Camin,Camax])
hold off
xlim([-Lmax Lmax])
ylim([-Lmax Lmax])
zlim([Camin Camax])
colormap(jet)
drawnow
end

% for p = 1:pf-1
% figure(pf)
% hold on
% plot(r,Ca(:,1,p))
% xlim([0 L])
% ylim([Cao Cab])
% end




