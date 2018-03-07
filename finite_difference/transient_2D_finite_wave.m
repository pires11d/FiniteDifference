clc
clf
clear all
close all

m = 8;
n = 8;
a = 1;
b = 1;
dx = 1/(m-1);
dy = 1/(n-1);
dt = 0.01;
tf = 10;

c = 10;

x = 0:dx:1;
y = 0:dy:1;
t = 0:dt:tf;

mf = length(x);
nf = length(y);
pf = length(t);

u = zeros(mf,nf,pf);

u(mf,:,:) = 0;      
u(:,nf,:) = 0;      
u(:,:,1) = 0;       

for p=1:pf-1
    for m = 2:mf-1
        for n = 2:nf-1
            X = (m-1)*dx;
            Y = (n-1)*dy;
            T = (p-1)*dt;
            phi = c*sin(a*pi*X)*sin(b*pi*Y)*cos(c*T);
            u(m,n,p+1) = u(m,n,p) + dt*(phi);
            u(m,n,p+2) = 2*u(m,n,p+1) - u(m,n,p) + (c^2)*dt*((u(m+1,n,p)-2*u(m,n,p)+u(m-1,n,p))/(dx^2)+(u(m,n+1,p)-2*u(m,n,p)+u(m,n-1,p))/(dy^2));
        end
    end
end

for p=1:pf-1
figure(1)
surf(x,y,u(:,:,p))
caxis([-1,1])
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 1])
ylim([0 1])
zlim([-1 1])
drawnow
end
