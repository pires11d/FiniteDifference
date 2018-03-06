clc
clf
clear all
close all

n = 9;
mf = n;
nf = n;
a = 1.0;
b = 1.0;
dx = a/(n-1);
dy = b/(n-1);
dt = 0.1;
tf = 30;

D = 0.008;

x = 0:dx:a;
y = 0:dy:b;
t = 0:dt:tf;

pf = length(t);

Ca = zeros(mf,nf,pf);
Cao = 100;
Cab = 0;
Camin = min(Cao,Cab);
Camax = max(Cao,Cab);

Ca(mf,:,:) = Cab;      
Ca(:,nf,:) = Cab;      
Ca(:,:,1) = Cao;       

for p=1:pf-1
    for m = 2:mf-1
        for n = 2:nf-1
            Ca(m,1,p+1)=Ca(m,2,p+1);
            Ca(1,n,p+1)=Ca(2,n,p+1);
            Ca(m,n,p+1)=Ca(m,n,p)+dt*D*((Ca(m+1,n,p)-2*Ca(m,n,p)+Ca(m-1,n,p))/(dx^2)+(Ca(m,n+1,p)-2*Ca(m,n,p)+Ca(m,n-1,p))/(dy^2));
        end
    end
    Ca(1,1,p+1)=(Ca(2,1,p+1)+Ca(1,2,p+1))/2;
end

for p=1:pf-1
figure(1)
surf(x,y,Ca(:,:,p))
hold on
surf(-x,y,Ca(:,:,p))
surf(x,-y,Ca(:,:,p))
surf(-x,-y,Ca(:,:,p))
caxis([Camin,Camax])
hold off
xlabel('X comprimento (m)')
ylabel('Y comprimento (m)')
zlabel('Concentracao (mol/L)')
xlim([-a a])
ylim([-b b])
zlim([Camin Camax])
colormap(jet)
drawnow
end






