% Nonisothermal PFR

clf
close all
clear all
clc

acetone_NIPFR_data

[X,Y] = ode45(@acetone_NIPFR_ODE,[0 Xf],[Vo,To]);

figure(1)
plot(X,Y(:,1),'b-')
grid on
ylabel('V (L)')
xlabel('X')
xlim([0 Xf])
ylim([0 5])

figure(2)
plot(X,Y(:,2),'r-')
grid on
ylabel('T (K)')
xlabel('X')
xlim([0 Xf])
ylim([800 1200])

