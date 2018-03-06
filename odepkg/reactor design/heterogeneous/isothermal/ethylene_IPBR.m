% Isothermal PBR

clf
close all
clear all
clc

ethylene_IPBR_data

[W,Y] = ode23(@ethylene_IPBR_ODE,[Wo Wf],[Xo yo]);

xx = Y(:,1);
yy = Y(:,2);
Pa = Pao*((thetaa-(a/a)*xx)./(1+e*xx)).*yy;
Pb = Pao*((thetab-(b/a)*xx)./(1+e*xx)).*yy;
Pc = Pao*((thetac+(c/a)*xx)./(1+e*xx)).*yy;
PI = Pao*(thetaI./(1+e*xx)).*yy;
f = (1+e*xx)./yy;
P = Po*yy;

figure(1)
plot(W,xx,'k-',W,yy,'m',W,f,'r')
grid on
ylabel('X')
xlabel('W (kg)')
xlim([0 Wf])
ylim([0 max(f)])
legend('X','y','f','Location','NorthWest')

figure(2)
plot(W,Pa,'b-',W,Pb,'g',W,Pc,'r',W,PI,'m',W,P,'k')
grid on
ylabel('Pi (atm)')
xlabel('W (kg)')
xlim([0 Wf])
ylim([0 Po])
legend('Pa','Pb','Pc','Pinert','P','Location','NorthEast')
