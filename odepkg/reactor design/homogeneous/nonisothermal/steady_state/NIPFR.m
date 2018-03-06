% Nonisothermal PFR

clf
close all
clear all
clc

NIPFR_data

[X,Y] = ode45(@NIPFR_ODE,[0:0.05:Xf],[Vo,To]);

figure(1)
plot(X,Y(:,1),'k-')
grid on
ylabel('V (L)')
xlabel('X')
xlim([0 Xf])
ylim([0 300])

figure(2)
plot(X,Y(:,2),'k-')
grid on
ylabel('T (K)')
xlabel('X')
xlim([0 Xf])
ylim([360 380])

Tf=Y(19,2)
Vf=Y(19,1)

