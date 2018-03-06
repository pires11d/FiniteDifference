% Isothermal PBR

clf
close all
clear all
clc

toluene_IPBR_data

[W,X] = ode23(@toluene_IPBR_ODE,[Wo Wf],Xo);

y = (1-a*W).^0.5;
Pa = Pao*(1-X).*y;
Pb = Pao*(thetab-X).*y;
Pc = Pao*(thetac+X).*y;

figure(1)
plot(W,X,'k-',W,y,'m')
grid on
ylabel('X')
xlabel('W (kg)')
xlim([0 Wf])
ylim([0 1])
legend('X','y','Location','NorthEast')

figure(2)
plot(W,Pa,'b-',W,Pb,'g',W,Pc,'r')
grid on
ylabel('Pi (atm)')
xlabel('W (kg)')
xlim([0 Wf])
ylim([0 20])
legend('Pa','Pb','Pc','Location','NorthEast')
