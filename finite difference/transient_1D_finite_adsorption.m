clc
clf
clear all
close all

m = 20;
L = 1;
dz = L/(m-1);
dt = 0.1;
tf = 100;

eps = 0.7;
rho = 1000;
u = 0.5;
K = 0.2;
D = 0.001;

z = 0:dz:L;
t = 0:dt:tf;

mf = length(z);
pf = length(t);

C = zeros(mf,pf);
Co = 0;
Ce = 1;

%% Condição Inicial
C(:,1) = Co;   
%% Condição de Contorno em z = 0
C(1,:) = Ce;      

%% DIFERENÇA FINITA
for p=1:pf-1
    for m = 2:mf-1
        %% Balanço de Massa
        adveccao_z = -u*(C(m,p)-C(m-1,p));
        difusao_z = (C(m+1,p)-2*C(m,p)+C(m-1,p))/(dz^2);
        C(m,p+1) = C(m,p) + (dt/(1+((1-eps)/eps)*rho*K))*(adveccao_z + difusao_z);
        %% Condição de Contorno em z = L
        C(mf,p+1) = (C(mf-1,p+1)+C(mf,p))/2;
    end
end

for p=1:pf-1
figure(1)
plot(z,C(:,p)/Ce)
title('Concentração Adsorvida Adimensional','FontSize',14)
xlabel('Z comprimento (m)')
ylabel('C/Co')
xlim([0 L])
ylim([0 Ce])
drawnow
end

figure(2)
mesh(t,z,C(:,:)/Ce)
title('Concentração Adsorvida Adimensional','FontSize',14)
xlabel('T tempo (s)')
ylabel('Z comprimento (m)')
zlabel('C/Co')
caxis([0,1])
view(-160,20)
hold off
xlim([0 tf])
ylim([0 L])
zlim([0 1])





