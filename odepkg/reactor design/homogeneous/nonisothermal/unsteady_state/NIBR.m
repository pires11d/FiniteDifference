% Nonisothermal Batch Reactor

clf
close all
clear all
clc

NIBR_data

[tt,Y] = ode45(@NIBR_ODE,[0:0.05:Xf],[to,To]);

figure(1)
plot(tt,Y(:,1),'k-')
grid on
ylabel('t (h)')
xlabel('X')
xlim([0 Xf])
ylim([0 3])

figure(2)
plot(tt,Y(:,2),'k-')
grid on
ylabel('T (K)')
xlabel('X')
xlim([0 Xf])
ylim([360 380])

