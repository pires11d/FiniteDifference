clf
close all
clear all
clc

% Nonisothermal PBR with pressure drop

NIPBR_data

%ODE:
[Z,Y] = ode23(@NIPBR_ODE,[0 zf],[Cao To Po Pao Xo Cbo Pbo Cco Pco]);

%Plots:
figure(1)
plot(Z,Y(:,2),'r-')
grid on
ylabel('T (K)')
xlabel('z (m)')
xlim([0 zf])
ylim([350 400])

figure(2)
plot(Z,Y(:,1),'b-',Z,Y(:,6),'g-',Z,Y(:,8),'m-')
grid on
ylabel('Ci (mol/L)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Cbo])
legend('Ca','Cb','Cc','Location','NorthEast')

figure(3)
plot(Z,Y(:,4),'b-',Z,Y(:,7),'g-',Z,Y(:,9),'r-',Z,Y(:,3),'k-')
grid on
ylabel('Pi (Pa)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Po])
legend('Pa','Pb','Pc','P','Location','NorthEast')

figure(4)
plot(Z,Y(:,5),'k-',Z,Y(:,3)/Po,'m')
grid on
ylabel('X')
xlabel('z (m)')
xlim([0 zf])
ylim([0 1])
legend('X','y','Location','SouthEast')
